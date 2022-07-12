%
% RunExamples for the main cases
addpaths;
clear all;
close all;

%% Set 2parameters
if ispc, SAVE_DIR = [getenv('USERPROFILE'), '\DataAnalyses\LearningDynamics']; else, SAVE_DIR = [getenv('HOME'), '/DataAnalyses/LearningDynamics']; end % Please keep this fixed, simply create a symlink ~/DataAnalyses pointing wherever you like                           
VERBOSE                         = 1;                                                                % indicator to print certain output
time_stamp                      = datestr(now, 30);
if ~exist('Params','var'), Params = [];     end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end

%% Load example definitions and let user select one example to run (choose OD or ODS)
Examples                        = LoadExampleDefinitions();
ExampleIdx                      = SelectExample(Params, Examples);

%% Get example parameters
Example                         = Examples{ExampleIdx};
sysInfo                        = Example.sysInfo;
solverInfo                     = Example.solverInfo;
obsInfo                        = Example.obsInfo;                                                 % move n to learn_info
obsInfo.VERBOSE                = VERBOSE;
obsInfo.SAVE_DIR               = SAVE_DIR;


% obs_info.obs_noise              = 0.0;
if obsInfo.obs_noise>0
  obsInfo.use_derivative     = true;
end

%% generate rhoT
obsInfo.MrhoT = 2000;

saveON= 0;
plotON = 0;
learnInfo.rhoLT = Generate_rhoT(sysInfo,obsInfo,solverInfo,saveON,plotON);% empirical pairwise distance

if obsInfo.obs_noise>0
  obsInfo.use_derivative     = true;
end

%% Perform learning and form the learnInfo
learnInfo.N = sysInfo.N;
learnInfo.d = sysInfo.d;
learnInfo.order = 2;
learnInfo.name = sysInfo.name;
learnInfo.jitter=1e-4;
learnInfo.Cov = 'Matern';
learnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation

learnInfo.v = 3/2; % 1/2, 3/2, 5/2, 7/2
learnInfo.hyp0 = log([1  1 NaN exp(1/2)]); %log([2 1]);% initialization of logsigma logomega, lognoise, m


if strcmp('ODS',learnInfo.name) 
   learnInfo.Ns = sysInfo.Ns;
   switch learnInfo.Ns
       case 1
           learnInfo.hyp0 = log([1  1 NaN 1/2 exp(1/2) exp(1/2)]); % initialization of logsigma logomega, lognoise, logk, m, Pa
       case 2
           learnInfo.hyp0 = log([1  1 NaN 1/2 exp(1/2) exp(1/2) exp(1/2)]); % initialization of logsigma logomega, lognoise, logk, m, Pa, Pb
       case 3    
           learnInfo.hyp0 = log([1  1 NaN 1/2 exp(1/2) exp(1/2) exp(1/2) exp(1/2)]); % initialization of logsigma logomega, lognoise, logk, m, Pa, Pb, Pc
   end
end

if obsInfo.obs_noise ~= 0
   learnInfo.hyp0(3) = log(1/2);
end



nT = 1;     %number of trials
errorphis = zeros(2,nT);      %store errors of phis in L-infinity and L2rhoT norms
hypparameters = zeros(length(learnInfo.hyp0),nT);   %store estimated hyperparameters


for k = 1:nT
    % Generate training data 
    [dxpath_test,xpath_test,dxpath_train, xpath_train]=Generate_training_data(sysInfo,obsInfo,solverInfo);

    dxpath_train= trajUnifNoiseAdditive(dxpath_train, obsInfo.obs_noise);

    saveON= 0;
    plotON = 0; 
    learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance 

    X_train=[];
    Y_train=[];

    %  For second order system:     X_train consists of [x_1,\cdots,x_N,v_1,\cdots,v_N]
    %                               Y_train consistss of [\dot v_1,..., \dot v_N]

    dt = (obsInfo.time_vec(2)-obsInfo.time_vec(1));
    timegap_train = floor(size(xpath_train,2)/3);


    for i = 1:obsInfo.M
        for j = 1:timegap_train:size(xpath_train,2)-1
            X_train = [X_train;xpath_train(:,j,i);dxpath_train(:,j,i);];

            accerl = (dxpath_train(:,j+1,i)-dxpath_train(:,j,i))./dt;
            Y_train = [Y_train;accerl];

        end
    end

    
    learnInfo.L= size(xpath_train,2)*size(xpath_train,3);
    learnInfo.X = X_train;
    learnInfo.Y = Y_train;

    learnInfo.xpath_train = xpath_train;
    learnInfo.dxpath_train = dxpath_train;
    
    learnInfo.hyp = learnInfo.hyp0;   
    
    Glik_hyp = @(hyp)Glik_MS(learnInfo,hyp);

    
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -600);
%     if flik > -10
%         [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -300);
%     end
    [~, ~,learnInfo] = Glik_MS(learnInfo,learnInfo.hyp);    
    

    
    hypparameters(:,k) = exp(learnInfo.hyp);
    if strcmp('OD',learnInfo.name) 
        hypparameters(4,k) = learnInfo.hyp(4);
    end
    if strcmp('ODS',learnInfo.name) 
        hypparameters(5:end,k) = learnInfo.hyp(5:end);
    end
    
    %% Make predictions on phi using these hyperparameters

    %Once the kernel is learned, store these once to avoid recomputing.
    learnInfo.invK = pinv(learnInfo.K);
    learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;

    fprintf('visualize the learning of the kernel...\n ');
    learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo);
    
    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
    
    [errorphis(1,k),errorphis(2,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range);
end


%% save files
filename=strcat(sysInfo.name,'M',num2str(obsInfo.M),'L',num2str(size(Y_train,1)/learnInfo.N/learnInfo.d/obsInfo.M),'sigma',num2str(obsInfo.obs_noise));
save(filename, 'sysInfo','obsInfo','learnInfo','hypparameters','errorphis');

%% Make predictions on trajectories using these hyperparameters
% if sysInfo.d ==1
%     visualize_trajs_1D(sysInfo,obsInfo,solverInfo,learnInfo);
%     
% elseif sysInfo.d ==2
%         visualize_trajs_2D(sysInfo,obsInfo,solverInfo,learnInfo);
% 
% end

% 
% 
% return




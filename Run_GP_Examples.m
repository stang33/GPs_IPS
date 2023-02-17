% (c) XXXX
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

%% Load example definitions and let user select one example to run
Examples                        = LoadExampleDefinitions();
ExampleIdx                      = SelectExample(Params, Examples);

%% Get example parameters
Example                         = Examples{ExampleIdx};
sysInfo                        = Example.sysInfo;
solverInfo                     = Example.solverInfo;
obsInfo                        = Example.obsInfo;                                                 % move n to learn_info
obsInfo.VERBOSE                = VERBOSE;
obsInfo.SAVE_DIR               = SAVE_DIR;


%% generate rhoT
obsInfo.MrhoT = 2000;

saveON= 0;
plotON = 0;
learnInfo.rhoLT = Generate_rhoT(sysInfo,obsInfo,solverInfo,saveON,plotON);% empirical pairwise distance

if obsInfo.obs_noise>0
  obsInfo.use_derivative     = true;
end


%% Compute error with multiple Learning Trials
learnInfo.N = sysInfo.N;
learnInfo.d = sysInfo.d;
learnInfo.order = sysInfo.ode_order;
learnInfo.name = sysInfo.name;
learnInfo.jitter=1e-4;
learnInfo.Cov = 'Matern';
learnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation
learnInfo.option = 'alldata';


learnInfo.v = 3/2; % 1/2, 3/2, 5/2, 7/2

if strcmp('ODS',learnInfo.name)   
   learnInfo.hyp0 = log([1  1 NaN 1/2 exp(1/2) exp(1/2) exp(1/2)]); % initialization of logsigma logomega, lognoise, logk, Pa, Pb, Pc
   if obsInfo.obs_noise ~= 0
       learnInfo.hyp0(3) = log(1/2);
   end
end

if strcmp('CSF',learnInfo.name) 
   learnInfo.hyp0 = log([1  1 NaN 1/2 3/2]); % initialization of logsigma logomega, lognoise, loga, logb
   if obsInfo.obs_noise ~= 0
       learnInfo.hyp0(3) = log(1/2);
   end
end

if strcmp('FM',learnInfo.name) 
   learnInfo.hyp0 = log([1  1 NaN 1 1]); % initialization of logsigma logomega, lognoise, loga, logb
   if obsInfo.obs_noise ~= 0
       learnInfo.hyp0(3) = log(1/2);
   end
end


nT = 10;     %number of trials
errorphis = zeros(2,nT);      %store errors of phis in L-infinity and L2rhoT norms
errortrajs_train = zeros(4,nT);     %store mean and std of trajectory error in training data
errortrajs_test = zeros(4,nT);      %store mean and std of trajectory error in testing data
hypparameters = zeros(length(learnInfo.hyp0),nT);   %store estimated hyperparameters

for k = 1:nT
    [dxpath_test,xpath_test,dxpath_train, xpath_train]=Generate_training_data(sysInfo,obsInfo,solverInfo);
    
    learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance 

    dxpath_train= trajUnifNoiseAdditive(dxpath_train, obsInfo.obs_noise);

    fprintf('\n Done! begin to learn ......');

    X_train=[];
    Y_train=[];

%  For second order system:     X_train consists of [x_1,\cdots,x_N,v_1,\cdots,v_N]
%                               Y_train consistss of [\dot v_1,..., \dot v_N]

    for i = 1:obsInfo.M
        for j = 1:size(xpath_train,2)
            X_train = [X_train;xpath_train(:,j,i)];
            if sysInfo.ode_order ==2
                Y_train = [Y_train;dxpath_train(sysInfo.d*sysInfo.N+1:end,j,i)];
            else
                Y_train = [Y_train;dxpath_train(:,j,i)];
            end
        end
    end
    
    learnInfo.L= size(xpath_train,2)*size(xpath_train,3);
    learnInfo.X = X_train;
    learnInfo.Y = Y_train;

    learnInfo.xpath_train = xpath_train;
    learnInfo.dxpath_train = dxpath_train;
    
    learnInfo.hyp = learnInfo.hyp0;
    Glik_hyp = @(hyp)Glik(learnInfo,hyp);

    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -50);
    [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);

    %Once the kernel is learned, store these once to avoid recomputing.
    learnInfo.invK = pinv(learnInfo.K);
    learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;
    
    hypparameters(:,k) = exp(learnInfo.hyp);
    if strcmp('ODS',learnInfo.name) 
        hypparameters(5:end,k) = learnInfo.hyp(5:end);
    end
    
    % %% visualize the kernel
    fprintf('visualize the learning of the kernel...\n ');

    learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo,'E');
    
    range = [0, learnInfo.rhoLT.edges(max(find(learnInfo.rhoLT.rdens~=0)))];
    
    [learnInfo,errorphis(1,k),errorphis(2,k)] = errornorms_phis(sysInfo,obsInfo,learnInfo,range,'E');
    result_train = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, learnInfo.xpath_train(:,1,:));
    errortrajs_train(:,k) = [result_train.train_traj_error result_train.prediction_traj_error]';
    result_test = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo,sysInfo.mu0());
    errortrajs_test(:,k) = [result_test.train_traj_error result_test.prediction_traj_error]';
end


%% save files
filename=strcat(sysInfo.name,'M',num2str(obsInfo.M),'L',num2str(length(obsInfo.time_vec)));
save(filename, 'sysInfo','obsInfo','learnInfo','errortrajs_train','errortrajs_test','hypparameters','errorphis');


%% visualize the' trajectory

if sysInfo.d ==1
    visualize_trajs_1D(sysInfo,obsInfo,solverInfo,learnInfo);
    
elseif sysInfo.d ==2
        visualize_trajs_2D(sysInfo,obsInfo,solverInfo,learnInfo);

end
 


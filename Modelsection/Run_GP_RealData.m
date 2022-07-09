% (c) XXXX
%% Test CSF with real data
addpaths;
clear all;
close all;

%% this data set contains frame 0 to frame 199. The orginnl dataset contains velocity
%% but we re-calculate the velocity after smoothing.

load 4601_Real_data.mat
% load 2201_Real_data.mat

% calculate the accerlation
d=2;
N=124;
dt = 1/3; % time gap


%% preprocessing the real data

% translatition to positive corordinates

x (1:2:end,:) = x(1:2:end,:)-min(min(x(1:2:end,:)));
x (2:2:end,:) = x(2:2:end,:)-min(min(x(2:2:end,:)));

% smooth data
for i=1:size(x,1)
    x(i,:)=smooth(x(i,:),10);
end




num2 = 1:2*N;

% normalize to [0,1]
x = (x)/max(max(x));

% compute the velocity
v = x(:,2:end)-x(:,1:end-1);
v = v./dt;


% extract position
x1 = x(num2,1:end-2);
% compute the accerlation
a1 = (v(num2,2:end)-v(num2,1:end-1))./dt;
%extract the velocity
v1 = v(num2,1:end-1);




sysInfo.name            = 'FM';                                                   % name of the dynamics
sysInfo.d               = 2;                                                                       % dimension for the state vector (in this case, opinion vector)
sysInfo.N               = N;                                                                      % # of agents
sysInfo.phi            = {@(r)FM_kernel(r,sysInfo.N), @(v,d) FM_ncforce(v,d)};                                               % energy based interaction
sysInfo.phi_type        = 'E';
sysInfo.K               = 1;                                                                       % # of types
sysInfo.ode_order       = 2;                                                                       % order of the ODE system
sysInfo.type_info       = ones(1, sysInfo.N);                                                     % function mapping agent index to its type index
sysInfo.RE              = [];   % energy based reulation on interactoin beween agent i and agent i'


% sysInfo.name            = 'CSF';                                                   % name of the dynamics
% sysInfo.d               = 2;                                                                       % dimension for the state vector (in this case, opinion vector)
% sysInfo.N               = N;                                                                      % # of agents
% %sysInfo.phi             = {@(r)CS_kernel(r), @(v,d) CS_force(v,d)};                                               % energy based interaction
% sysInfo.phi_type        = 'v';
% sysInfo.K               = 1;                                                                       % # of types
% sysInfo.ode_order       = 2;                                                                       % order of the ODE system
% sysInfo.type_info       = ones(1, sysInfo.N);                                                     % function mapping agent index to its type index
% sysInfo.RE              = [];                                                                      % energy based reulation on interactoin beween agent i and agent i'



sysInfo.domain          = [floor(min(min(x1))), ceil(max(max(x1)))];
sysInfo.type = 1;
sysInfo.type_info = ones(sysInfo.N,1);

obsInfo.M               = 1;                                                                      % # trajectories with random initial conditions for learning interaction kernel
obsInfo.rho_T_histedges    = linspace(0,sysInfo.domain(2)-sysInfo.domain(1),1000);  % a rather arbitrary set of bins to track estimators of \rho^L_T




%%% chose the simulation time

training_opt = 1; % set as 0 if no training

tc = 40; % timescale
test_time_vec =0:tc/197:tc;


% state vector Y=[X,V]
xpath = [x1;v1];
% vector Z =[V,\dot V]
dxpath = [v1;a1];

time = [1 29];
xpath_train = xpath(:,time);
dxpath_train = dxpath(:,time);

saveON= 0;
plotON = 1;
learnInfo.rho_emp = rho_empirical(xpath_train,sysInfo,obsInfo,saveON,plotON);% empirical pairwise distance

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

learnInfo.N = sysInfo.N;
learnInfo.d = sysInfo.d;
learnInfo.order = sysInfo.ode_order;
learnInfo.name = sysInfo.name;
learnInfo.L= size(xpath_train,2)*size(xpath_train,3);
learnInfo.X = X_train;
learnInfo.Y = Y_train;
learnInfo.jitter=1e-5;
learnInfo.Cov = 'Matern';
learnInfo.dtraj ='false'; % do not compute dtraj during trajectory computation

learnInfo.v = 3/2; % 1/2, 3/2, 5/2, 7/2

if strcmp('CSF',learnInfo.name) 
   % learnInfo.hyp0 = log([1  1 0.001 1 1]); % initialization of logsigma logomega, lognoise, loga, logb
   learnInfo.hyp0 = log([rand(1,4) 0.001 rand(1,2)]);
end

if strcmp('FM',learnInfo.name) 
   learnInfo.hyp0 = log([rand(1,4) 0.001 rand(1,2)]);
   % learnInfo.hyp0 = log([rand(1,4) NaN rand(1,2)]);
end


learnInfo.xpath_train = xpath_train;
learnInfo.dxpath_train = dxpath_train;


%% train hyperameters
learnInfo.hyp = learnInfo.hyp0;   

learnInfo.option = 'subset';  % doesn't work for ODS
learnInfo.Nsub = 2;
learnInfo.sub = randsample(1:learnInfo.N,learnInfo.Nsub);

Glik_hyp = @(hyp)Glik(learnInfo,hyp);

if training_opt
    [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -100);
end

% [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -600);
% [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
% learnInfo.CoefM = pinv(learnInfo.K)*learnInfo.Ym; % compute the coeficient matrix

learnInfo.option = 'alldata';
[~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);

%Once the kernel is learned, store these once to avoid recomputing.
learnInfo.invK = pinv(learnInfo.K);
learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;

learnInfo = visualize_phis_realdata(sysInfo,obsInfo,learnInfo,'E');
learnInfo = visualize_phis_realdata(sysInfo,obsInfo,learnInfo,'A');

%%

traj_hat = construct_traj(learnInfo,test_time_vec,xpath_train(:,1));
if strcmp('subset',learnInfo.option)
    if training.opt
        filename=strcat(sysInfo.name,'Realdatalearning100','S',num2str(learnInfo.Nsub),'.mat');
    else
        filename=strcat(sysInfo.name,'Realdatalearning','S',num2str(learnInfo.Nsub),'.mat');
    end
else
    if training.opt
        filename=strcat(sysInfo.name,'Realdatalearning100','.mat');
    else
        filename=strcat(sysInfo.name,'Realdatalearning','.mat');
    end
end
save(filename, 'sysInfo','obsInfo','learnInfo','traj_hat');




%%visualize the trajectory
visualize_trajs_2D_realdata(sysInfo, xpath(:,1:96),traj_hat);
%% Test FM with real data
addpaths;
clear all;
close all;

load 4601_Real_data.mat
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

sysInfo.domain          = [floor(min(min(x1))), ceil(max(max(x1)))];
sysInfo.type = 1;
sysInfo.type_info = ones(sysInfo.N,1);

obsInfo.M               = 1;                                                                      % # trajectories with random initial conditions for learning interaction kernel
obsInfo.rho_T_histedges    = linspace(0,sysInfo.domain(2)-sysInfo.domain(1),1000);  % a rather arbitrary set of bins to track estimators of \rho^L_T

sysInfo.flagxi          = 0;

%%% chose the simulation time

training_opt = 1; % set as 0 if no training

if training_opt
%   test_time_vec = 0:0.2:19;
    tc = 38; % timescale
    test_time_vec =0:0.4:38;
else
    test_time_vec =0:3/197:3;
end


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



% %% train hyperameters
% learnInfo.hyp = learnInfo.hyp0;   
% 
% learnInfo.option = 'subset';  % doesn't work for ODS
% learnInfo.Nsub = 12;
% learnInfo.sub = randsample(1:learnInfo.N,learnInfo.Nsub);
% 
% Glik_hyp = @(hyp)Glik(learnInfo,hyp);
% 
% if training_opt
%     [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -100);
% end
% 
% % [learnInfo.hyp,flik,i] = minimize(learnInfo.hyp, Glik_hyp, -600);
% % [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
% % learnInfo.CoefM = pinv(learnInfo.K)*learnInfo.Ym; % compute the coeficient matrix
% 
% % learnInfo.option = 'alldata';
% [~, ~,learnInfo] = Glik(learnInfo,learnInfo.hyp);
% 
% %Once the kernel is learned, store these once to avoid recomputing.
% learnInfo.invK = pinv(learnInfo.K);
% learnInfo.invKprodYm = learnInfo.invK * learnInfo.Ym;
% 
% learnInfo = visualize_phis_realdata(sysInfo,obsInfo,learnInfo,'E');
% learnInfo = visualize_phis_realdata(sysInfo,obsInfo,learnInfo,'A');
% 
% %%
% 
% traj_hat = construct_traj(learnInfo,test_time_vec,xpath_train(:,1));
% if strcmp('subset',learnInfo.option)
%     if training.opt
%         filename=strcat(sysInfo.name,'Realdatalearning100','S',num2str(learnInfo.Nsub),'.mat');
%     else
%         filename=strcat(sysInfo.name,'Realdatalearning','S',num2str(learnInfo.Nsub),'.mat');
%     end
% else
%     if training.opt
%         filename=strcat(sysInfo.name,'Realdatalearning100','.mat');
%     else
%         filename=strcat(sysInfo.name,'Realdatalearning','.mat');
%     end
% end
% save(filename, 'sysInfo','obsInfo','learnInfo','traj_hat');
% 
% 
% 
% 
% %%visualize the trajectory
% visualize_trajs_2D_realdata(sysInfo, xpath(:,1:96),traj_hat);
%% using Sindy to learn dynamics
%  Matlab SINDy package is available at faculty.washington.edu/sbrunton/sparsedynamics.zip.

% if do comparison
Comparison_with_SINDy= true;

if Comparison_with_SINDy
	home_path = [pwd filesep];
    addpath(genpath([home_path '/SINDy_utils/']));
    
    % prepare for the X and dot X
    
    Sindyfit_tic =tic;
    
    dxdata = dxpath_train;
    xdata  = xpath_train;
%     dxdata     = xdata(:,2:end,:) - xdata(:,1:end-1,:);
%     dt = obs_info.time_vec(2)-obs_info.time_vec(1);
%     train_dxdata_1 = reshape(dxdata./dt,size(dxdata,1),[])'; % Data of dot X, size LM x Nd
    
    train_xdata_1 = reshape(xdata,size(xdata,1),[])'; % Data of X, size LM x Nd
    train_dxdata_1 = reshape(dxdata,size(dxdata,1),[])'; % Data of dot X, size LM x Nd
    

    %% pool Data  (i.e., build library of nonlinear time series)
    
    polyorder = 2;
    usesine   = 1;
    Theta = poolData(train_xdata_1,size(train_xdata_1,2),polyorder,usesine); % the defaul dictionary is mutivariable polys with sines and cosines
    m = size(Theta,2);% size of dictionary
    
    %% compute Sparse regression: sequential least squares
    % lambda is our sparsification knob. in our example, sparsity is not obvious,
    %  have to set lambda very small, large lambda leads to zero solution
    lambda = 0.00005;
    Xi = sparsifyDynamics(Theta,train_dxdata_1,lambda,size(train_xdata_1,2)); % coefficient matrix, size m x Nd
    fprintf('\n--- The elapse time for Running Sindy is: %10.4e',toc(Sindyfit_tic));
    
    
    
    %% use sparse regression coefficient as RHS
    
    SindyPred_tic =tic;
    RHS = @(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine);
    
    
    %% prepare data for X and dot X

	dxdata = dxpath_train;
    xdata  = xpath_train;
%     dt = obsInfo.time_vec(2)-obsInfo.time_vec(1);
%     dxdata     = xdata(:,2:end,:) - xdata(:,1:end-1,:);
%     train_dxdata = reshape(dxdata./dt,size(dxdata,1),[]); % size: Nd x LM

    train_xdata = reshape(xdata,size(xdata,1),[]); % Data of X, size LM x Nd
    train_dxdata = reshape(dxdata,size(dxdata,1),[]); % Data of dot X, size LM x Nd
    
    dxdata_true = train_dxdata;
    obsfut_true = reshape(xpath(:,1:96),size(dxdata,1),[]);
    dxdata_fut = reshape(dxpath(:,1:96),size(dxdata,1),[]);
    
    for l=1: size(train_xdata,2)
%         dxdata_true(:,l)                   = RHSfn_2nd_ncf(t,train_xdata(:,l),N,sysInfo.phi{1},sysInfo.phi{2},sysInfo.phi_type);                                                  % for future usage when the sys_info (especially type_info) is updated during the time integration
        hatdxdata_true(:,l)                = predict_ode_interp(train_xdata(:,l),learnInfo);                                               % for future usage when the sys_info (especially type_info) is updated during the time integration 
    end

    for l=1: size(obsfut_true,2)
%         dxdata_fut(:,l)                   = eval_rhs( obsfut_true(:,l), sys_info);                                                    % for future usage when the sys_info (especially type_info) is updated during the time integration
        hatdxdata_fut(:,l)                   = predict_ode_interp(obsfut_true(:,l),learnInfo);                                                    % for future usage when the sys_info (especially type_info) is updated during the time integration
    end
%
    fprintf('\n--- RMSE of our algorithm of fitting fphi for Training ICs on training time interval is: %10.4e',norm(hatdxdata_true-dxdata_true)./norm(dxdata_true));
    fprintf('\n--- RMSE of our algorithm of fitting fphi for Training ICS on future time interval is: %10.4e',norm(hatdxdata_fut-dxdata_fut)./norm(dxdata_fut));



    %% Training error for Training ICs
    SINDy_dxdata = train_xdata;
    SINDy_dxdata_fut = obsfut_true;
    for j=1:size(train_xdata,2)
        SINDy_dxdata(:,j)= RHS(0,train_xdata(:,j));
    end
    
    for j=1:size(obsfut_true,2)
        SINDy_dxdata_fut(:,j)= RHS(0,obsfut_true(:,j));
    end
    
    perf_sindy_train = norm(SINDy_dxdata(sysInfo.d*sysInfo.N+1:end,:)-train_dxdata(sysInfo.d*sysInfo.N+1:end,:),2)./norm(train_dxdata(sysInfo.d*sysInfo.N+1:end,:));
    perf_sindy_train_fut = norm(SINDy_dxdata_fut(sysInfo.d*sysInfo.N+1:end,:)-dxdata_fut(sysInfo.d*sysInfo.N+1:end,:))./norm(dxdata_fut(sysInfo.d*sysInfo.N+1:end,:));
    fprintf('\n--- RMSE error of fitting fphi using Sindy on training time interval is: %10.4e',perf_sindy_train);
    fprintf('\n--- RMSE error of fitting fphi using Sindy on future time interval is: %10.4e',perf_sindy_train_fut);
    
    
    %% Parameters for Trajectory prediction
    sysInfo.T_f = 38;
    tspan = [0,sysInfo.T_f];
    options = odeset('RelTol',1e-5,'AbsTol',1e-6);
        
    %% Trajectory prediction for training ICs on both training time interval and future time interval
    ICs = learnInfo.xpath_train(:,1,1);
%     for i=1:size(ICs,2)
%         sol_SINDy = ode15s(RHS,tspan,ICs(:,i),options);
%         xpath_SINDy_train(:,:,i) = deval(sol_SINDy,obsInfo_Ltest.time_vec);
%         xpath_SINDy_test(:,:,i) = deval (sol_SINDy,obsInfo_Ltest_fut.time_vec);
%     end

    time_vec = 0:0.4:38;    
    
    sol_SINDy = ode15s(RHS,tspan,ICs(:,1),options);
	xpath_SINDy(:,:,1) = deval(sol_SINDy,time_vec);
    
    sup_err_sindy               = traj_norm(obsfut_true(1:sysInfo.d*sysInfo.N,:),     xpath_SINDy(1:sysInfo.d*sysInfo.N,:,1),    'Time-Maxed', sysInfo);
    

    fprintf('\n------------------- Overall time for computing trajecotry errors: %.2f',toc(SindyPred_tic));
%     
%     
    fprintf('\n--- SINDy: trajectory prediciton is: %10.4e', sup_err_sindy );
    
%% visualize the results for one trajectory prediction

	visualize_trajs_2D_realdata(sysInfo, xpath(:,1:96),xpath_SINDy);

    %% Train a neural_network to approximate  RHS and use NN to predict the dynamics
    %
    NNfit_tic=tic;

	dxdata = dxpath_train;
    xdata  = xpath_train;
%     dt = obsInfo.time_vec(2)-obsInfo.time_vec(1);
%     dxdata     = xdata(:,2:end,:) - xdata(:,1:end-1,:);
%     train_dxdata = reshape(dxdata./dt,size(dxdata,1),[]); % Data of dot X, % Nd x LM

    train_xdata = reshape(xdata,size(xdata,1),[]); % Data of X, %Nd x LM
    train_dxdata = reshape(dxdata,size(dxdata,1),[]); % Data of dot X, % Nd x LM
    
    dxdata_true = train_dxdata;
    
    obsfut_true = reshape(xpath(:,1:96),size(xdata,1),[]);
    dxdata_fut = reshape(dxpath(:,1:96),size(dxdata,1),[]);

     
    %% construct and train neural net
    %net = fitnet(200); % two hidden layer of 25 neurons each
    net = feedforwardnet([40,20]);
    net = train(net,train_xdata,train_dxdata);
    
    %% evaluate the performance of fitting fphi on the training data set
    net_dxdata= net(train_xdata);
    net_dxdata_fut = net(obsfut_true);
    perf_train = norm(net_dxdata(sysInfo.d*sysInfo.N+1:end,:)-train_dxdata(sysInfo.d*sysInfo.N+1:end,:))/norm(train_dxdata(sysInfo.d*sysInfo.N+1:end,:));
    perf_train_fut = norm(net_dxdata_fut(sysInfo.d*sysInfo.N+1:end,:)-dxdata_fut(sysInfo.d*sysInfo.N+1:end,:))/norm(dxdata_fut(sysInfo.d*sysInfo.N+1:end,:));
    fprintf('\n--- The elapse time for fitting neural network is: %10.4e',toc(NNfit_tic));
    fprintf('\n--- RMSE error of fitting fphi using FNN fitting fphi on training time interval is: %10.4e',perf_train);
    fprintf('\n--- RMSE error of fitting fphi using FNN fitting fphi on future time interval is: %10.4e',perf_train_fut);
    
    
    %% use neural network as RHS
    
    NNPred_tic =tic;
    RHS = @(t,x)RHS_NN(net,x);
    tspan = [0,sysInfo.T_f];
    options = odeset('RelTol',1e-5,'AbsTol',1e-6);

    time_vec = 0:0.4:38;    
    
    sol_NN = ode15s(RHS,tspan,ICs(:,1),options);
	xpath_NN(:,:,1) = deval(sol_NN,time_vec);

    sup_err               = traj_norm(obsfut_true(1:sysInfo.d*sysInfo.N,:),     xpath_NN(1:sysInfo.d*sysInfo.N,:,1),    'Time-Maxed', sysInfo);


	fprintf('\n------------------- Overall time for computing trajecotry errors: %.2f',toc(NNPred_tic));
    fprintf('\n--- FNN: trajectory prediciton is: %10.4e', sup_err);

	visualize_trajs_2D_realdata(sysInfo,xpath(:,1:96),xpath_NN);
   
    save(strcat(sysInfo.name,'_SINDy_NN_err'),'sup_err_sindy','sup_err');
    fprintf('\ndone.\n');
    

    %% compare the performance 
    [~,M_t] = group_polarisation(xpath(:,1:96),sysInfo);
    [~,M_gp] = group_polarisation(traj_hat(:,1:96),sysInfo);
    [~,M_SINDy] = group_polarisation(xpath_SINDy(:,1:96),sysInfo);
    [~,M_NN] = group_polarisation(xpath_NN(:,1:96),sysInfo);

    f1 = figure;
    tiledlayout(2,2)
    
    ax1 = nexttile;
    plot(ax1,1:96,M_t)
    title(ax1,'True Traj')
    xlim(ax1,[1 96])
    xlabel(ax1,'t')
    ylabel(ax1,'|M(t)|')
    
    ax2 = nexttile;
    plot(ax2,1:96,M_gp)
    title(ax2,'GP Approx')
    xlim(ax2,[1 96])
    xlabel(ax2,'t')
    ylabel(ax2,'|M(t)|')
    
    ax3 = nexttile;
    plot(ax3,1:96,M_SINDy)
    title(ax3,'SINDy Approx')
    xlim(ax3,[1 96])
    xlabel(ax3,'t')
    ylabel(ax3,'|M(t)|')
    
    ax4 = nexttile;
    plot(ax4,1:96,M_NN)
    title(ax4,'FNN Approx')
    xlim(ax4,[1 96])
    xlabel(ax4,'t')
    ylabel(ax4,'|M(t)|')
    
    f2 = figure;
    tiledlayout(2,2)
    
    ax5 = nexttile;
    h1 = histogram(M_t,'BinWidth',0.01);
    title(ax5,'True Traj')
    xlabel(ax5,'|M(t)|')
    xlim(ax5,[0.05 0.3])
    
    
    ax6 = nexttile;
    h2 = histogram(M_gp,'BinWidth',0.01);
    title(ax6,'GP Approx')
    xlabel(ax6,'|M(t)|')
    xlim(ax6,[0.05 0.3])

    ax7 = nexttile;
    h3 = histogram(M_SINDy,'BinWidth',0.01);
    title(ax7,'SINDy Approx')
    xlabel(ax7,'|M(t)|')
    xlim(ax7,[0.05 0.3])
    
    ax8 = nexttile;
    h4 = histogram(M_NN,'BinWidth',0.01);
    title(ax8,'FNN Approx')
    xlabel(ax8,'|M(t)|')
    xlim(ax8,[0.05 0.3])
end


%% this code is to visualize the trajectories in the stochastic dynamic system
% kerEst
% rhoLT
% rhoLT_empirical

% Authored by Jinchao Feng and Sui Tang

function visualize_training_trajs_1D(sysInfo,obsInfo,solverInfo, learnInfo)

%% set up the plot information


 
    
% changing plot_info
plot_info.scrsz                 = [1, 1, 1920, 1080];                                        
plot_info.legend_font_size      = 33;
plot_info.legend_font_name      = 'Helvetica';
plot_info.colorbar_font_size    = 33;
plot_info.title_font_size       = 44;
plot_info.title_font_name       = 'Helvetica';
plot_info.axis_font_size        = 44;
plot_info.axis_font_name        = 'Helvetica';
plot_info.tick_font_size        = 39;
plot_info.tick_font_name        = 'Helvetica';
plot_info.traj_line_width       = 1.0;
plot_info.phi_line_width        = 1.5;
plot_info.phihat_line_width     = 1.5;
plot_info.rhotscalingdownfactor = 1;
plot_info.showplottitles        = false;
plot_info.display_phihat        = false;
plot_info.display_interpolant   = true;
plot_info.T_L_marker_size       = 40;
plot_info.line_styles           = {'-' '-.', '--', ':'};                                        % traj. line styles
plot_info.type_colors           = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};
plot_info.marker_style          = {'s', 'd', 'p', 'h', 'x', '+', 'v', '^', '<', '>'};            
plot_info.marker_size           = plot_info.T_L_marker_size;                                                   
plot_info.marker_edge_color     = {'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'};
plot_info.marker_face_color     = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};
plot_info.plot_name = sysInfo.name;
plot_info.make_movie =0;
if plot_info.make_movie, plot_info.movie_name = sprintf('%s/%s_movie_%s', sys_info.name, num2str(sysInfo.T),num2str(obsInfo.M)); end

% go through the Initial Conditions from the training data set
num_trajs = 1;
trajs                = cell(1, num_trajs); 


%% for initial data drawn from traning data
data = learnInfo.xpath_train;
IC=data(:,1,1);

N         = sysInfo.N;         % number of agents   
d         = sysInfo.d; 
order = sysInfo.ode_order;
dN = N*d*learnInfo.order;
myODE =  @(t,x) RHSfn(t,x,N,sysInfo.phi{1});

if order ==1&&strcmp(learnInfo.name,'ODNS')
    if sysInfo.Noption == 1
        myODE =  @(t,x) RHSfn_NS1(t,x,N,sysInfo.phi{1});
    else
        myODE =  @(t,x) RHSfn_NS(t,x,N,sysInfo.phi{1});
    end
end

if order ==2&&strcmp(learnInfo.name,'CSNS')
    if sysInfo.Noption == 1
        myODE =  @(t,y) RHSfn_2nd_NS1(t,y,N,sysInfo.phi{1},sysInfo.phi_type);
    else
        myODE =  @(t,y) RHSfn_2nd_NS(t,y,N,sysInfo.phi{1},sysInfo.phi_type);
    end
end

train_time_vec=obsInfo.time_vec;
test_time_vec = [train_time_vec(1:end-1) linspace(train_time_vec(end),sysInfo.T_f,50)];% prediction on [0, T_f]                                                                     % final time the system will reach steady state
L = length(test_time_vec); %

traj_true            = zeros (dN,L,size(IC, 2)); % true traj

m = 1;
  
sol_true = ode45(myODE,solverInfo.time_span,IC(:,m),solverInfo.option); % solu from adaptive solver
traj_true(:,:,m)  =deval(sol_true,test_time_vec) ;% Nd x steps

trajs{1}      = traj_true(1:d*N,:,:);


%% visualize  1D trajectories

% prepare the window size
if isfield(plot_info, 'scrsz') && ~isempty(plot_info.scrsz), scrsz = plot_info.scrsz; else, scrsz = get(groot,'ScreenSize'); end
% prepare the figure window
traj_fig                       = figure('Name', 'Traj (1D): training', 'NumberTitle', 'off', 'Position', ...
  [scrsz(3)*1/8 + scrsz(3) * 6/48, scrsz(4)*1/8, scrsz(3)*3/4, scrsz(4)*3/4]);
% prepare the y range
y_min                          = zeros(1, length(trajs));
y_max                          = zeros(1, length(trajs));
for ind = 1 : length(trajs)
  traj                         = trajs{ind};
  y_min(ind)                   = min(min(traj));
  y_max(ind)                   = max(max(traj));
end
y_min                          = min(y_min);
y_max                          = max(y_max);
if length(trajs) == 1, handleAxes = gobjects(1); else, handleAxes = gobjects(length(trajs)/2, 2); end
y_range                        = y_max - y_min;
y_min                          = y_min - 0.1 * y_range;
y_max                          = y_max + 0.1 * y_range;
vline_y                        = linspace(y_min, y_max, size(train_time_vec,2));
vline_x                        = train_time_vec(end) * ones(size(vline_y));
x_min                          = min(train_time_vec);
x_max                          = max(test_time_vec);
for ind = 1 : length(trajs)
  if length(trajs) == 1, sp_handle = subplot(1, 1, ind); else,  sp_handle = subplot(length(trajs)/2, 2, ind); end
  traj                         = trajs{ind};

  
  %% plot 1D trajectory
  
   for k = 1 : sysInfo.type
      index     = sysInfo.type_info == k;
      agents_traj = traj(index, :);
      plot(sp_handle, test_time_vec, agents_traj, 'LineWidth', plot_info.traj_line_width, 'Color', plot_info.type_colors{k}, ...
          'LineStyle', plot_info.line_styles{1});
      if k == 1, hold on; end
   end


  plot(vline_x, vline_y, '-.k');
  l_handle                     = plot(NaN, NaN, ['k' plot_info.line_styles{1}]);  % dummy lines for legend
  hold off;
  axis([x_min, x_max, y_min, y_max]);
  sp_handle.FontSize           = plot_info.tick_font_size;
  sp_handle.FontName           = plot_info.tick_font_name;   
  if mod(ind, 2) == 1
    ylabel('Coord. $1$', 'Interpreter', 'latex', 'FontSize', plot_info.axis_font_size, 'FontName', plot_info.axis_font_name);
  end
  if ind == length(trajs) - 1 || ind == length(trajs)
    xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', plot_info.axis_font_size, 'FontName', plot_info.axis_font_name);
 end
 leg_handle                 = legend(l_handle, {'$\mathbf{x}_i(t)$'});
 set(leg_handle, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', plot_info.legend_font_size);
  row_ind                      = floor((ind - 1)/2) + 1;
  col_ind                      = mod(ind - 1, 2) + 1;   

  handleAxes(row_ind, col_ind) = sp_handle;
end

saveas(traj_fig, [plot_info.plot_name '_traj'], 'fig'); 



end

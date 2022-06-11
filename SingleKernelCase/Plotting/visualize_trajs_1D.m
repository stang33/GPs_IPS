%% this code is to visualize the trajectories in the stochastic dynamic system
% kerEst
% rhoLT
% rhoLT_empirical

% (c) XXXX

function visualize_trajs_1D(sysInfo,obsInfo,solverInfo, learnInfo)

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
num_trajs = 4;
trajs                = cell(1, num_trajs); 


%% for initial data drawn from traning data
data = learnInfo.xpath_train;
IC=data(:,1,1);
result               = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, IC);
trajs{1}             = result.traj_true;
trajs{2}             = result.traj_hat;

%% for initial data drawn from test data
IC =  sysInfo.mu0();
result               = construct_and_compute_traj(sysInfo,obsInfo,solverInfo,learnInfo, IC);
trajs{3}             = result.traj_true;
trajs{4}             = result.traj_hat;
train_time_vec             = result.train_time_vec;
prediction_time_vec        = result.prediction_time_vec;
time_vec=[train_time_vec prediction_time_vec];



%% visualize  1D trajectories

% prepare the window size
if isfield(plot_info, 'scrsz') && ~isempty(plot_info.scrsz), scrsz = plot_info.scrsz; else, scrsz = get(groot,'ScreenSize'); end
% prepare the figure window
traj_fig                       = figure('Name', 'Traj (1D): True Vs. Learned', 'NumberTitle', 'off', 'Position', ...
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
x_max                          = max(prediction_time_vec);
for ind = 1 : length(trajs)
  if length(trajs) == 1, sp_handle = subplot(1, 1, ind); else,  sp_handle = subplot(length(trajs)/2, 2, ind); end
  traj                         = trajs{ind};

  
  %% plot 1D trajectory
  
   for k = 1 : sysInfo.type
      index     = sysInfo.type_info == k;
      agents_traj = traj(index, :);
      plot(sp_handle, time_vec, agents_traj, 'LineWidth', plot_info.traj_line_width, 'Color', plot_info.type_colors{k}, ...
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
  if ind == 1 || ind == 3
    leg_handle                 = legend(l_handle, {'$\mathbf{x}_i(t)$','$\hat\mathbf{x}_i(t)$'});
  elseif ind == 2 || ind == 4
    leg_handle                 = legend(l_handle, {'$\mathbf{x}_i(t)$','$\hat\mathbf{x}_i(t)$'});

  elseif ind == 5 || ind == 6
    leg_handle                 = legend(l_handle, {'$\mathbf{x}_i^{LN}(t)$'});
  end
 set(leg_handle, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', plot_info.legend_font_size);
  row_ind                      = floor((ind - 1)/2) + 1;
  col_ind                      = mod(ind - 1, 2) + 1;   

  handleAxes(row_ind, col_ind) = sp_handle;
end

saveas(traj_fig, [plot_info.plot_name '_traj'], 'fig'); 

%Generate movie



end







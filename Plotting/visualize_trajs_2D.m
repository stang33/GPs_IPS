%% this code is to visualize the trajectories in the stochastic dynamic system
% kerEst
% rhoLT
% rhoLT_empirical

% (c) XXXX

function visualize_trajs_2D(sysInfo,obsInfo,solverInfo, learnInfo)


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
plot_info.traj_line_width       = 2.0;
plot_info.phi_line_width        = 1.5;
plot_info.phihat_line_width     = 1.5;
plot_info.rhotscalingdownfactor = 1;
plot_info.showplottitles        = false;
plot_info.display_phihat        = false;
plot_info.display_interpolant   = true;
plot_info.arrow_thickness          = 1.5;                                                           % thickness of the arrow body
plot_info.arrow_head_size          = 0.8;                                                           % size of the arrow head
plot_info.arrow_scale              = 0.05;  
plot_info.T_L_marker_size       = 10;

plot_info.line_styles           = {'-', '-.', '--', ':'};                                        % traj. line styles
plot_info.type_colors           = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};
plot_info.marker_style          = {'s', 'd', 'p', 'h', 'x', '+', 'v', '^', '<', '>'};
plot_info.marker_size           = plot_info.T_L_marker_size;
plot_info.marker_edge_color     = {'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'};
plot_info.marker_face_color     = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};
plot_info.plot_name = sysInfo.name;
plot_info.make_movie = 0;


% go through the Initial Conditions from the training data set
num_trajs = 4;
trajs                = cell(1, num_trajs);
dtrajs               = cell(1, num_trajs);


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
time_vec=[train_time_vec(1:end) prediction_time_vec];


%% visualize  2D trajectories

% prepare the window size
if isfield(plot_info, 'scrsz') && ~isempty(plot_info.scrsz), scrsz = plot_info.scrsz; else, scrsz = get(groot,'ScreenSize'); end
% prepare the figure window
traj_fig                       = figure('Name', 'Traj (2D): True Vs. Learned', 'NumberTitle', 'off', 'Position', ...
    [scrsz(3)*1/8 + scrsz(3) * 6/48, scrsz(4)*1/8, scrsz(3)*3/4, scrsz(4)*3/4]);

color_output                   = construct_color_items(sysInfo.type, train_time_vec(end), time_vec);
cmap                           = color_output.cmap;
plot_info.c_vecs               = color_output.c_vecs;
clabels                        = color_output.clabels;
cticks                         = color_output.cticks;



% split the trajectories

X_coords1                        = cell(length(trajs), sysInfo.d);

the_mins                       = zeros(length(trajs), sysInfo.d);
the_maxs                       = zeros(length(trajs), sysInfo.d);

for ind = 1 : length(trajs)
    traj1                         = trajs{ind}; % true traj: Nd x timesteps

    for d_ind = 1 : sysInfo.d
        X_cd1                       = traj1(d_ind : sysInfo.d : end - (sysInfo.d - d_ind), :); %  dth -component of trajectory

        the_mins(ind, d_ind)       = min(min(X_cd1(:)));%,min(X_cd2(:)));
        the_maxs(ind, d_ind)       = max(max(X_cd1(:)));%,max(X_cd2(:)));
        X_coords1{ind, d_ind}       = X_cd1;

    end
end
x_min                          = min(the_mins(:, 1));
x_max                          = max(the_maxs(:, 1));
y_min                          = min(the_mins(:, 2));
y_max                          = max(the_maxs(:, 2));


if length(trajs) == 1, handleAxes = gobjects(1); else, handleAxes = gobjects(1,length(trajs)/2); end
T_loc                          = find(time_vec == train_time_vec(end));

for ind = 1 : length(trajs)
    if length(trajs) == 1, sp_handle = subplot(1, 1, ind); else, sp_handle = subplot(length(trajs)/2, 2, ind); end

    X_cds1 = X_coords1(ind, :);

    %% plot 2D trajectory
    set(traj_fig, 'CurrentAxes', sp_handle);
    line_style_ind=1;
    for k = 1 : sysInfo.type
        agents_Ck                  = find(sysInfo.type_info == k);
        N_k                        = length(agents_Ck);
        for agent_ind = 1 : N_k
            agent                    = agents_Ck(agent_ind);
            if line_style_ind < 3
                c1_at_t_1                = X_cds1{1}(agent, :);

                c2_at_t_1                = X_cds1{2}(agent, :);

            elseif line_style_ind == 3
                c1_at_t_1                = [X_cds1{1}(agent, :), NaN * ones(1, plot_info.c_len)];

                c2_at_t_1               = [X_cds1{2}(agent, :), NaN * ones(1, plot_info.c_len)];

            end
      
           p_handle                 = patch([c1_at_t_1, NaN], [c2_at_t_1, NaN], [plot_info.c_vecs{k}, NaN], 'EdgeColor', ...
                    'interp', 'LineStyle', plot_info.line_styles{line_style_ind}, 'LineWidth', plot_info.traj_line_width);
            
                
      % display the uncertainty region covariance

            if strcmp(sysInfo.name, 'PredatorPrey1stOrder') || strcmp(sysInfo.name, 'PredatorPrey2ndOrder') || strcmp(sysInfo.name, 'PredatorPrey1stOrderSplines')
                if k == 1, set(p_handle, 'EdgeAlpha', 0.25); end
            end
            if k == 1 && agent_ind == 1, hold on; end
        end
    end
    
    
    
    
    plot(X_coords1{ind, 1}(:, T_loc), X_coords1{ind, 2}(:, T_loc), 'o', 'MarkerSize', plot_info.T_L_marker_size, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    l_handle                   = plot(NaN, NaN, ['k' plot_info.line_styles{1}]);                    % dummy lines for legend
    
    hold off;
    if ind == 1
        xticks                     = get(sp_handle, 'XTick');
        delta                      = xticks(2) - xticks(1);
    end
    xtickformat('%.1f'); ytickformat('%.1f');
    if sysInfo.d == 2, axis([x_min, x_max + delta, y_min, y_max]); else, axis([x_min, x_max, y_min, y_max, z_min, z_max]); end
    sp_handle.FontSize           = plot_info.tick_font_size;
    sp_handle.FontName           = plot_info.tick_font_name;
    if length(trajs) == 1 || mod(ind, 2) == 0
        colormap(cmap);
        cbh                        = colorbar('YTickLabel', clabels, 'YTick', cticks, 'Location', 'East');
        set(cbh, 'FontSize', plot_info.colorbar_font_size);
    end
    if mod(ind, 2) == 1
        if sysInfo.d == 2
            ylabel('Coord. $2$', 'FontSize', plot_info.axis_font_size, 'Interpreter', 'latex', 'FontName', plot_info.axis_font_name);
        else
            zlabel('Coord. $3$', 'FontSize', plot_info.axis_font_size, 'Interpreter', 'latex', 'FontName', plot_info.axis_font_name);
        end
    end
    if ind == length(trajs) - 1 || ind == length(trajs)
        if sysInfo.d == 2
            xlabel('Coord. $1$', 'FontSize', plot_info.axis_font_size, 'Interpreter', 'latex', 'FontName', plot_info.axis_font_name);
        else
            xlabel('Coord. $1$', 'FontSize', plot_info.axis_font_size, 'Interpreter', 'latex', 'FontName', plot_info.axis_font_name);
            ylabel('Coord. $2$', 'FontSize', plot_info.axis_font_size, 'Interpreter', 'latex', 'FontName', plot_info.axis_font_name);
        end
    end
    row_ind                      = floor((ind - 1)/2) + 1;
    col_ind                      = mod(ind - 1, 2) + 1;
    if ind == 1 || ind == 3
        leg_handle                 = legend(l_handle, {'$\mathbf{x}_i(t)$'});
    elseif ind == 2 || ind == 4
        leg_handle                 = legend(l_handle, {'$\hat\mathbf{x}_i(t)$'});
    elseif ind == 5 || ind == 6
        leg_handle                 = legend(l_handle, {'$\mathbf{x}_i^{LN}(t)$'});
    end
    set(leg_handle, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', plot_info.legend_font_size);
    handleAxes(row_ind, col_ind) = sp_handle;
end

saveas(traj_fig, [plot_info.plot_name '_traj'], 'fig');

if plot_info.make_movie, plot_info.movie_name = sprintf('%s/%s/%s/%s_movie',pwd, 'outputs',sysInfo.name,sysInfo.name); end

if plot_info.make_movie 
 make_traj_multiD_animation(trajs, dtrajs, time_vec, sysInfo, plot_info)
end



end

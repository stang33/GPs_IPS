function learnInfo = visualize_phis(sysInfo,obsInfo,learnInfo)
% function visualize_phis(learningOutput, sys_info, obs_info, plot_info)

% (c) M. Zhong, S.Tang(JHU)




%% set up the plot information
plot_info.scrsz                 = [1, 1, 1920, 1080];                                        
plot_info.legend_font_size      = 33;
plot_info.legend_font_name      = 'Helvetica';
plot_info.colorbar_font_size    = 33;
plot_info.title_font_size       = 44;
plot_info.title_font_name       = 'Helvetica';
plot_info.axis_font_size        = 30;
plot_info.axis_font_name        = 'Helvetica';
plot_info.tick_font_size        = 30;
plot_info.tick_font_name        = 'Helvetica';
plot_info.traj_line_width       = 1.0;
plot_info.phi_line_width        = 1.5;
plot_info.phihat_line_width     = 1.5;
plot_info.rhotscalingdownfactor = 1;
plot_info.showplottitles        = false;
plot_info.display_phihat        = false;
plot_info.display_interpolant   = true;
plot_info.T_L_marker_size       = 40;
plot_info.line_styles           = {'-', '-.', '--', ':'};                                        % traj. line styles
plot_info.type_colors           = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};
plot_info.marker_style          = {'s', 'd', 'p', 'h', 'x', '+', 'v', '^', '<', '>'};            
plot_info.marker_size           = plot_info.T_L_marker_size;                                                   
plot_info.marker_edge_color     = {'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'};
plot_info.marker_face_color     = {'b', 'r', 'c', 'g', 'm', 'y', 'k', 'k', 'k', 'k'};


% screen size issue
if isfield(plot_info, 'scrsz'), scrsz = plot_info.scrsz; else, scrsz = get(groot,'ScreenSize'); end
% do the energy based phi's first

% prepare the window to plot the true and learned interactions  
  phi_fig = figure('Name', 'PhiEs: True Vs. Learned', 'NumberTitle', 'off', 'Position', ...
            [scrsz(3)*1/8, scrsz(4)*1/8, scrsz(3)*3/4, scrsz(4)*3/4]);
  
%% load the parameters
rho_emp = learnInfo.rho_emp;
plot_info.plot_name = strcat(sysInfo.name,'_results_phis');
one_knot                =  rho_emp.edges;                 % find out the knot vector
range     = [one_knot(1), one_knot(200)];                                 % plot the uniform learning result first; use the knot vectors as the range
r           = one_knot(1):(one_knot(2)-one_knot(1)):2*range(2);                                % refine this knot vector, so that each sub-interval generated a knot has at least 7 interior points, so 3 levels of refinement
%r(r<range(1) | r>range(2))  = [];
%r                           = r(r>=0 & r<=rho_emp.edgesSupp(end) );

%% for phi in the state variable 
phi                       = sysInfo.phi{1}(r);
[phi_mean,phi_cov]=phi_pos(r,learnInfo);
cov=diag(phi_cov);
cov(cov<0)=0;
phi_cov_diag = sqrt(cov);
y_min               = min([min(phi-2*phi_cov_diag),min(phi_mean-2*phi_cov_diag)])*1;                                                                       % find out the range for y values, f1 might have blow up, just use f_hat
y_max               = max([max(phi+2*phi_cov_diag),max(phi_mean+2*phi_cov_diag)])*1;
if y_max<y_min+10*eps, y_max=y_min+1; y_min=y_min-1; end
set(groot, 'CurrentFigure', phi_fig);                                               % plot everything on the true_vs_learned window

edges               = obsInfo.rho_T_histedges;  % Estimated \rho's
edges_idxs_fine          = find(range(1) <= edges & edges<range(2));
%centers_fine             = (edges(edges_idxs_fine(1):edges_idxs_fine(end)-1) + edges(edges_idxs_fine(1)+1:edges_idxs_fine(end)))/2;
% downsampling of hist data by 50
edges_idxs =edges_idxs_fine(1:5:end);
%centers =centers_fine(1:20:end);
centers =edges(edges_idxs);
histdata1            = rho_emp.rdens(edges_idxs(1:end));                    % this is the "true" \rhoLT from many MC simulations

%% save the result posterior mean
learnInfo.r = r;
learnInfo.phi = phi_mean;
learnInfo.phi_cov = phi_cov_diag;
      
%% plot the density of rho
yyaxis right                                                                                % display \rho^L_T and its estimator
axesHandle1         = gca();
%histHandle1         = plot(centers,histdata1,'k-','LineWidth',1);    hold on;
hist1       = bar(centers,histdata1,'FaceColor', [.8 .8 .8]);    hold on;
hist1.FaceAlpha = 0.5;
%hist1               = fill(centers(([1 1:end end])),[0 histdata1 0], ...
%    'k','EdgeColor','none','FaceAlpha',0.2);

axesHandle1.FontSize           = plot_info.tick_font_size;
axesHandle1.FontName           = plot_info.tick_font_name; 
axis([range(1), range(2), 0, 2]);
axis tight;

%ylabelHandle = ylabel(axesHandle1,'$\hat\rho^L_T, \rho^{L,M}_T$','Interpreter','latex','FontSize',AXIS_FONT_SIZE);
%set(ylabelHandle,'VerticalAlignment','bottom');
%set(ylabelHandle,'Position',get(ylabelHandle,'Position').*[0.92,0.33,1]);





%% plot the kernel
yyaxis left
display interaction kernels
plotHandle1=plot(r, phi, '-k', 'LineWidth', plot_info.phi_line_width);hold on
plotHandle2=plot(r, phi_mean,  '-b', 'LineWidth', plot_info.phihat_line_width);               % plot the learned phi


%% display the uncertainty region covariance 
rconf = [r r(end:-1:1)] ;
yconf = [phi_mean'+2*phi_cov_diag' phi_mean(end:-1:1)'-2*phi_cov_diag(end:-1:1)'];
p = fill(rconf,yconf,'red');
p.FaceColor = [0.7 0.8 0.9];      
p.EdgeColor = 'none';   
set(p,'facealpha',.5)

    
legendHandle2=legend([plotHandle1,plotHandle2,hist1], ...
        {'$\phi$','$\hat{\phi}_{mean}$','$\hat\rho^L_T$'});
    %ylabel(['$\phi, \hat{\phi}, \hat{\phi}_{reg}$'],'Interpreter','latex','FontSize',AXIS_FONT_SIZE);
    set(legendHandle2,'Location', 'NorthEast','Interpreter','latex','FontSize',plot_info.legend_font_size);



% for opinion dynamics axis([range(1), range(2), -0.6, 1.5]); 
axis([range(1), range(2), y_min, y_max]); 


% set up a uniform x-range, tighten the y-range
xlabel('r (pairwise distances)','FontSize',plot_info.axis_font_size);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + 2*ti(1);
bottom = outerpos(2) + 1.5*ti(2);
ax_width = outerpos(3) - 2*ti(1) - 2*ti(3);
ax_height = outerpos(4) - 1.7*ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

mkdir('./figures');
saveas(phi_fig, [plot_info.plot_name '_phi'], 'fig');






end
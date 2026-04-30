% Authored by Jinchao Feng and Sui Tang
% make_figure.m
% Publication-quality single-panel trajectory-error figures for the JSC revision.
% Produces 4 PDFs in ../images/:
%   ODNSd1_traj_error.pdf   -- 1st-order opinion dynamics, d=1
%   ODNSd2_traj_error.pdf   -- 1st-order opinion dynamics, d=2
%   CSNS1_traj_error.pdf    -- Cucker-Smale, cut-off kernel  (paper Example 1)
%   CSNS2_traj_error.pdf    -- Cucker-Smale, smooth kernel   (paper Example 2)
% Each shows relative max-time predict-window error vs M, one curve per noise
% level, log y-axis with +/- 1 std shaded band.
clear; close all;

M_list  = 1:10;
data_dir = '../trajectory_data';
imdir    = '../images';
if ~exist(imdir,'dir'), mkdir(imdir); end

% Each entry:
%   .file_prefix = on-disk .mat prefix (matches what's in trajectory_data/)
%   .out_name    = output PDF base name (matches paper's image naming)
%   .title       = title shown in panel
%   .nsr         = list of noise-to-signal ratios for this system
%   .dim_token   = 'd1' / 'd2' for the dimension in filenames
panels = { ...
  struct('file_prefix','ODNS1','dim_token','d1','out_name','ODNSd1_traj_error', ...
         'nsr',[0 0.25 0.5 0.75 1], ...
         'title','ODNS, $d=1$ (step kernel)'); ...
  struct('file_prefix','ODNS2','dim_token','d2','out_name','ODNSd2_traj_error', ...
         'nsr',[0 0.25 0.5 0.75 1], ...
         'title','ODNS, $d=2$ (step kernel)'); ...
  struct('file_prefix','CSNS2','dim_token','d2','out_name','CSNS1_traj_error', ...
         'nsr',[0 0.1 0.2 0.3 0.4 0.5], ...
         'title','CSNS1: cut-off kernel $\phi(r)=\mathbf{1}_{[0,0.5]}(r)$'); ...
  struct('file_prefix','CSNS1','dim_token','d2','out_name','CSNS2_traj_error', ...
         'nsr',[0 0.025 0.05 0.075 0.1], ...
         'title','CSNS2: smooth kernel $\phi(r)=1/(1+r^2)^{1/4}$') };

floor_val = 1e-6;
linestyles = {'-','--','-.',':','-','--','-.'};
markers    = {'o','s','^','d','v','>','<'};

for p = 1:numel(panels)
    P = panels{p};
    nsrs = P.nsr;
    cmap = viridis_like(numel(nsrs));

    % --- load data ---
    mu = NaN(numel(M_list), numel(nsrs));
    sd = NaN(numel(M_list), numel(nsrs));
    for mi = 1:numel(M_list)
        for ki = 1:numel(nsrs)
            fname = sprintf('%s/%sN100%sM%dL6sigma%s.mat', data_dir, ...
                            P.file_prefix, P.dim_token, M_list(mi), num2str(nsrs(ki)));
            if exist(fname,'file')
                S = load(fname);
                vals = S.errortrajs_test(3,:);  % predict-window error per trial (fresh-IC)
                mu(mi,ki) = mean(vals);
                sd(mi,ki) = std (vals);
            else
                fprintf('WARN: missing %s\n', fname);
            end
        end
    end

    % --- figure ---
    fig = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
    ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax,'YScale','log','FontSize',10,'TickLabelInterpreter','latex', ...
           'GridAlpha',0.2,'MinorGridAlpha',0.07,'Layer','top');
    h = gobjects(1,numel(nsrs));
    for ki = 1:numel(nsrs)
        m = mu(:,ki); s = sd(:,ki);
        upper = m + s; lower = max(m - s, floor_val);
        patch(ax,'XData',[M_list, fliplr(M_list)], ...
                 'YData',[upper', fliplr(lower')], ...
                 'FaceColor',cmap(ki,:),'FaceAlpha',0.07,'EdgeColor','none', ...
                 'HandleVisibility','off');
        h(ki) = plot(ax, M_list, m, ...
                     'LineStyle',linestyles{ki},'Marker',markers{ki}, ...
                     'Color',cmap(ki,:),'MarkerFaceColor',cmap(ki,:), ...
                     'MarkerSize',4.5,'LineWidth',1.5, ...
                     'DisplayName',sprintf('$\\sigma_{\\mathrm{nsr}}=%g$',nsrs(ki)));
    end
    xlabel(ax,'$M$','Interpreter','latex','FontSize',11);
    ylabel(ax,'$\mathrm{err}_{\mathrm{pred}}$','Interpreter','latex','FontSize',11);
    title(ax, P.title, 'Interpreter','latex','FontSize',11);
    xlim(ax,[0.5 10.5]); xticks(ax,1:10);
    lg = legend(ax, h,'Interpreter','latex','Location','best', ...
                'FontSize',8.5,'EdgeColor',[0.7 0.7 0.7],'NumColumns',1);
    lg.ItemTokenSize = [16 9];

    out_pdf = fullfile(imdir, [P.out_name '.pdf']);
    out_png = fullfile(imdir, [P.out_name '.png']);
    exportgraphics(fig, out_pdf, 'ContentType','vector');
    exportgraphics(fig, out_png, 'Resolution',220);
    close(fig);
    fprintf('saved %s\n', out_pdf);
end

%% --- helper: viridis-like sequential colormap ---
function C = viridis_like(n)
    base = [
        0.267 0.005 0.329
        0.282 0.140 0.458
        0.254 0.265 0.530
        0.207 0.372 0.553
        0.164 0.471 0.558
        0.128 0.567 0.551
        0.135 0.659 0.518
        0.267 0.749 0.441
        0.478 0.821 0.318
        0.741 0.873 0.150
        0.993 0.906 0.144 ];
    idx = round(linspace(1, size(base,1), n));
    C = base(idx, :);
end

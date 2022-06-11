function line_handles = plot_one_point_multiD(fig_handle, plot_handle, X_cds, dX_cds, time_ind, arrow_len, is_at_T, sys_info, plot_info)
% function plot_one_point_multiD(fig_handle, plot_handle, X_cds, dX_cds, time_ind, arrow_len, is_at_T, sys_info, plot_info)

% (c) XXXX

set(fig_handle, 'CurrentAxes', plot_handle);
line_handles                 = gobjects(1, sys_info.type);
for k = 1 : sys_info.type
  agents_Ck                  = find(sys_info.type_info == k);
  N_k                        = length(agents_Ck);
  for agent_ind = 1 : N_k
    agent                    = agents_Ck(agent_ind);
    if ~is_at_T
      if sys_info.d == 2
        ph = plot(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), plot_info.marker_style{k}, 'MarkerSize', plot_info.marker_size, ...
          'MarkerEdgeColor', plot_info.marker_edge_color{k}, 'MarkerFaceColor', plot_info.marker_face_color{k});
      else
        ph = plot3(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), X_cds{3}(agent, time_ind), plot_info.marker_style{k}, 'MarkerSize', ...
          plot_info.marker_size, 'MarkerEdgeColor', plot_info.marker_edge_color{k}, 'MarkerFaceColor', plot_info.marker_face_color{k});        
      end
    else
      if sys_info.d == 2
        ph = plot(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), 'o', 'MarkerSize', plot_info.marker_size, ...
          'MarkerEdgeColor', plot_info.marker_edge_color{k}, 'MarkerFaceColor', plot_info.marker_face_color{k});
      else
        ph = plot3(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), X_cds{3}(agent, time_ind), 'o', 'MarkerSize', ...
          plot_info.marker_size, 'MarkerEdgeColor', plot_info.marker_edge_color{k}, 'MarkerFaceColor', plot_info.marker_face_color{k});        
      end
    end
    if agent_ind == 1, line_handles(k) = ph; end
    if k == 1 && agent_ind == 1, hold on; end
    if sys_info.d == 2
      quiver(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), arrow_len(agent) * dX_cds{1}(agent, time_ind), arrow_len(agent) * dX_cds{2}(agent, time_ind), ...
        'LineWidth', plot_info.arrow_thickness, 'Color', plot_info.marker_face_color{k}, 'ShowArrowHead', 'on', 'MaxHeadSize', plot_info.arrow_head_size, ...
        'AutoScale', 'off');
    else
      quiver3(X_cds{1}(agent, time_ind), X_cds{2}(agent, time_ind), X_cds{3}(agent, time_ind), arrow_len(agent) * dX_cds{1}(agent, time_ind), ...
        arrow_len(agent) * dX_cds{2}(agent, time_ind), arrow_len(agent) * dX_cds{3}(agent, time_ind), 'LineWidth', plot_info.arrow_thickness, ...
        'Color', plot_info.marker_face_color{k}, 'ShowArrowHead', 'on', 'MaxHeadSize', plot_info.arrow_head_size, 'AutoScale', 'off');      
    end
  end
end
end
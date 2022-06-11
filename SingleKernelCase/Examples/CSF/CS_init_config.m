function y_init = CS_init_config(L_x, L_y,d, N, kind)
%
% function y_init = LJ_truncated_init_config(L_x, d, N, kind)
%

% (c) XXXX

% generate the initial configuration based on kind
switch kind
    case 1
        y_init_position = uniform_dist(d, N, 'rectangle', [-L_x, L_x]);
        y_init_velocity = uniform_dist(d, N, 'rectangle', [-L_y, L_y]);
        y_init = [y_init_position(:);y_init_velocity(:)];
   
end

return
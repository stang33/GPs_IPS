function y_init = FM_init_config(L_x,R_x,L_v,R_v, d, N, kind)
%
% function y_init = LJ_truncated_init_config(L_x, d, N, kind)
%

% (c) XXXX

% generate the initial configuration based on kind
switch kind
    case 1
        y_init_position = uniform_dist(d, N, 'rectangle', [L_x, R_x]);
        y_init_velocity = uniform_dist(d, N, 'rectangle', [L_v, R_v]);
        y_init = [y_init_position(:);y_init_velocity(:)];
   
end

return
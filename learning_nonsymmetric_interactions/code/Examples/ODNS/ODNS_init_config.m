function y_init = ODNS_init_config(d, N, kind)

%

% Authored by Jinchao Feng and Sui Tang


switch kind
  case 1
    if d ~= 1
      error('OD_init_config:exception', 'For 1D example, d has to be 1!!');
    end
    y_init = uniform_dist(d, N, 'line', [0, 10]);
    y_init = y_init(:);

  case 2
    if d == 1
      error('OD_init_config:exception', 'For 2D example, d has to be 2!!');
    end 
    y_init = uniform_dist(d, N, 'rectangle', [0, 4]);
    y_init = y_init(:);
  otherwise
end
end

 


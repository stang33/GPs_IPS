function out = getuij(data, i, j, l, m, d)

    p_i = [data(d * i - 1, l, m), data(d * i, l, m)];
    p_j = [data(d * j - 1, l, m), data(d * j, l, m)];

    out = p_i - p_j;


end
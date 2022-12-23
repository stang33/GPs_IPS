function out = getSpatialCoord(data, i, j, m, l)

    d = 2; %ambient dimension of particles
    out = data((d * i - 1) + (j - 1), l, m);

end
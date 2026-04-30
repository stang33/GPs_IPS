function val = phi_piecewise_const(r, coef_L, rmax)

    nbasis = length(coef_L);
    dr = rmax / nbasis;

    % Preserve original shape of r
    sz = size(r);

    % Convert r to bin indices
    idx = floor(r(:) ./ dr) + 1;
    idx = min(max(idx, 1), nbasis);

    % Lookup values
    val = coef_L(idx);

    % Restore shape
    val = reshape(val, sz);

end
function rhoVect = GetRhoVect(n, M, L, omega, delta, ds)

    tildeN = (n * (n - 1) * M * L) / 2;

    %Store the needed rhos for the sparse inverse.
    rhoVect = zeros(tildeN - 1,1);
    
    for i = 1 : tildeN - 1
        rhoVect(i,1) = KForFullRExp(ds(i+1), ds(i), omega, delta, n);
    end


end
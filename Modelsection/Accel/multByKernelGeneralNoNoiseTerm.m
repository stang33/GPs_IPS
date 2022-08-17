function zTilde = multByKernelGeneralNoNoiseTerm(Z, sigma, n, D, M, L, U_re, P_r, P_c, rhoVect, newRawIndices, delta, ls, ms)

%General N,M,L.
rhoVectTemp = (n^2 / delta^2) * rhoVect;
capN = n * D * M * L;
tildeN = (n * (n - 1) * M * L) / 2;

%Calculate g1 given decomp.
g1 = zeros(tildeN, 1);

for k = 1 : tildeN
    
    l = ls(k);
    m = ms(k);

    sumTotal = 0;
    cTerm = @(j_k) (D-1) * (m-1) * n * L + (D - 1) * (l - 1) * n + n * j_k;
    
    for j_k = 0 : D - 1
        sumTotal = sumTotal + U_re(P_r(k,1) + cTerm(j_k), P_r(k,2) - (m-1)*n*L - (l-1)*n)   *   (Z(2 * P_r(k,1) - 1 + j_k,1) - Z(2 * P_r(k,2) - 1 + j_k,1));
    end

    g1(k,1) = sumTotal;


end


%Calculate g2 given decomp.
%This depends on the kernel used.
g2 = zeros(tildeN,1);

%Edge case 1: tildeN
g2(tildeN) = g1(tildeN) * sqrt(1 - rhoVectTemp(tildeN - 1)^2) / (n/delta);

for k = tildeN - 1 : -1 : 2    
    g2(k) = sqrt(1 - rhoVectTemp(k - 1)^2) * g1(k) / (n/delta) + (1 / sqrt(1 - rhoVectTemp(k)^2)) * (rhoVectTemp(k) * g2(k + 1) * sqrt(1 - rhoVectTemp(k-1)^2));
end

%Edge case 2: 1
g2(1) = rhoVectTemp(1) / sqrt(1 - rhoVectTemp(1)^2) * g2(2) + g1(1) / (n/delta);



%Calculate g3 given decomp.
%This depends on the kernel used.    
g3 = zeros(tildeN,1);    
g3(1) = 1/(n/delta) * g2(1);

for k = 2 : tildeN    
    g3(k) = sqrt(1 - rhoVectTemp(k-1)^2) / (n/delta) * g2(k) + rhoVectTemp(k-1) * g3(k-1);

end




%Now calculate g4 as the sparse product.
sparseProd = zeros(capN,1);    

for m = 1 : M
    for l = 1 : L
        for i = 1 : n
            for j = 1 : D
    
                row = U_re((j-1)*n + i + (l-1)*n*D + (m-1)*n*D*L, :);

                %Delete 0s from same particle dists.
                row(row==0) = [];
        
                %Retreive correct indices from decomp.
                indexSet = P_c(i+ (l-1)*n + (m-1)*L*n,:);    
                rawIndexSet = newRawIndices(i+ (l-1)*n+ (m-1)*L*n,:);
        
                sum = 0;
                for k = 1 : n - 1
                    sum = sum + (g3(indexSet(k)) * row(rawIndexSet(k)));
                end
        
                %D*(i-1)+j is the k entry of the tildeN vector.
                sparseProd(D*(i-1)+j+ (l-1)*n*D+ (m-1)*n*D*L) = sum;
        
            end
        end
    end
end

%Output is ready.
zTilde = sparseProd;

end

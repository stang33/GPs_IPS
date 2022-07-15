function zTilde = multByKernel(Z, sigma, n, D, M, L, U_re, P_r, P_c, rhoVect, newRawIndices)

    capN = n * D * M * L;
    tildeN = (n * (n - 1) * M * L) / 2;

    %Calculate g1 given decomp.
    g1 = zeros(tildeN, 1);

    for k = 1 : tildeN
    
        m = 1;
        l = 1;
        cTerm = @(j_k) (D-1) * (m-1) * n * L + (D - 1) * (l - 1) * n + n * j_k;
    
        sumTotal = 0;
        for j_k = 0 : D - 1
            sumTotal = sumTotal + U_re(P_r(k,1) + cTerm(j_k), P_r(k,2) - (m-1)*n*L - (l-1)*n)   *   (Z(2 * P_r(k,1) - 1 + j_k,1) - Z(2 * P_r(k,2) - 1 + j_k,1));
        end
    
        g1(k,1) = sumTotal;
    
    end

    %Calculate g2 given decomp.
    %This depends on the kernel used.
    g2 = zeros(tildeN,1);
    
    %Edge case 1: tildeN
    g2(tildeN) = g1(tildeN) * sqrt(1 - rhoVect(tildeN - 1)^2);
    
    for k = tildeN - 1 : -1 : 2    
        g2(k) = sqrt(1 - rhoVect(k - 1)^2) * g1(k) + (1 / sqrt(1 - rhoVect(k)^2)) * (rhoVect(k) * g2(k + 1) * sqrt(1 - rhoVect(k-1)^2));
    end
    
    %Edge case 2: 1
    g2(1) = rhoVect(1) / sqrt(1 - rhoVect(1)^2) * g2(2) + g1(1);
    
    

    %Calculate g3 given decomp.
    %This depends on the kernel used.    
    g3 = zeros(tildeN,1);    
    g3(1) = g2(1);
    
    for k = 2 : tildeN    
        g3(k) = sqrt(1 - rhoVect(k-1)^2) * g2(k) + rhoVect(k-1) * g3(k-1);
    
    end
    
    %Now calculate g4 as the sparse product.
    sparseProd = zeros(capN,1);    
    
    for i = 1 : n
        for j = 1 : D

            row = U_re((j-1)*n + i, :);
    
            %Delete 0s from same particle dists.
            row(row==0) = [];
    
            %Retreive correct indices from decomp.
            indexSet = P_c(i,:);    
            rawIndexSet = newRawIndices(i,:);
    
            sum = 0;
            for k = 1 : n - 1
                sum = sum + (g3(indexSet(k)) * row(rawIndexSet(k)));
            end
    
            %D*(i-1)+j is the k entry of the tildeN vector.
            sparseProd(D*(i-1)+j) = sum;
    
        end
    end

    
    %Output is ready.
    zTilde = sparseProd + sigma * Z;
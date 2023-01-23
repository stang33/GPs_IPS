function K = ConstructKernelMatrix(distanceData, outerData, omega, delta, n, D, M, L, nu, partialOmega)

    capN = n * D * M * L;
    tildeN = (n * (n - 1) * M * L) / 2;
    
    U_s = zeros(capN,tildeN);
    
    %Build U_s with entries that are the distance pairs.
    colNum = 1;
    for m = 1 : M
        for l = 1 : L
            for i = 1 : n-1
                for j = i + 1 : n
            
                    U_s(((D * i - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*i + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(outerData,i,j,l,m,D);
                    U_s(((D * j - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*j + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(outerData,j,i,l,m,D);
            
                    colNum = colNum + 1;
                    
                end
            end
        end
    end

    %Now, we make R_s.
    
    %First, get distances.
    
    ds = zeros(tildeN,1);
    entry = 1;
    for m = 1 : M
        for l = 1 : L
            for i = 1 : n
                for j = i + 1 : n
                    ds(entry) = norm(getuij(distanceData,i,j,l,m,D));
                    entry = entry + 1;
                end
            end
        end
    end
    
    
    [ds, indices] = sort(ds);
    

    %Swap columns.
    U_s = U_s(:,indices);
    R_s = zeros(tildeN, tildeN);
    for i = 1:tildeN
        for j = 1:tildeN  

            R_s(i,j) = MaternKernel(ds(i),ds(j), omega, delta,n, nu, partialOmega);

        end
    end

    %Return result.
    K = U_s * R_s * U_s';
    
 

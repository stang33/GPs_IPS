function [U_s, R_s, UnsortU, sortingInd] = TotalDecompForDebug(data, omega, delta, n, D, M, L, learnInfo)

    capN = n * D * M * L;
    tildeN = (n * (n - 1) * M * L) / 2;
    
    U_s = zeros(capN,tildeN);
    
    
    %Try to build U_s with entries that are the distance pairs
    colNum = 1;
    for m = 1 : M
        for l = 1 : L
            for i = 1 : n-1
                for j = i + 1 : n
            
                    U_s(((D * i - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*i + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(data,i,j,l,m);
                    U_s(((D * j - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*j + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(data,j,i,l,m);
            
                    colNum = colNum + 1;
                    
                end
            end
        end
    end

    %Now, we try to make R_s.
    
    %First, make the distances.
    
    ds = zeros(tildeN,1);
    entry = 1;
    for m = 1 : M
        for l = 1 : L
            for i = 1 : n
                for j = i + 1 : n
                    ds(entry) = norm(getuij(data,i,j,l,m));
                    entry = entry + 1;
                end
            end
        end
    end
    
    
    [ds, indices] = sort(ds);
    

    %Swap cols
    %2 lines debug
    UnsortU = U_s;
    sortingInd = indices;
    U_s = U_s(:,indices);
    
    
    if learnInfo.v == 1/2

        R_s = zeros(tildeN, tildeN);
        
        for i = 1:tildeN
            for j = 1:tildeN
        
                R_s(i,j) = KForFullRExp(ds(i),ds(j), omega, delta,n);
        
            end
        end

    elseif learnInfo.v == 5/2

    
        gamma = 1 / omega;
        sigma = delta;
        %TODO: What are these params? The above is maybe right?
        lambda = sqrt(5) / gamma;
        q = 16/3 * sigma^2 * lambda^5;
    
    
        Wis = zeros(3,3,tildeN);
        Wis(1,1,1) = sigma^2;
        Wis(1,2,1) = 1;
        Wis(1,3,1) = -sigma^2*lambda^2/3;
        Wis(2,2,1) = sigma^2*lambda^2/3;
        Wis(2,3,1) = 1;
        Wis(3,1,1) =  Wis(1,3,1);
        Wis(3,3,1) = sigma^2*lambda^4;
    
        Gis = zeros(3,3,tildeN);
    
        for i = 1 : tildeN
        
            W = zeros(3,3);
            di = ds(i);
            W(1,1) = (exp(-2*lambda*di) * (3 + 6 * lambda * di + 6 * lambda^2 * di^2 + 4 * lambda^3 * di^3 + 2 * lambda^4 * di^4) - 3)/ (-4 * lambda^5);
            W(1,2) = (exp(-2*lambda*di) * di^4) / 2;
            W(1,3) = (exp(-2*lambda*di) * (1 + 2*lambda*di + 2*lambda^2*di^2 + 4*lambda^3*di^3-2*lambda^4*di^4) - 1)/(4*lambda^3);
            W(2,1) = W(1,2);
            W(3,1) = W(1,3);
            W(2,2) = (exp(-2*lambda*di) * (1 + 2 * lambda * di + 2 * lambda^2 * di^2 - 4 * lambda^3 * di^3 + 2 * lambda^4 * di^4) - 1)/(-4*lambda^3);
            W(2,3) = (exp(-2*lambda*di) * di^2 * (4 -4*lambda*di+lambda^2*di^2))/ 2;
            W(3,2) = W(2,3);
            W(3,3) = (exp(-2*lambda*di) * (-3 + 10 * lambda^2 * di^2 - 22 * lambda^2 * di^2 + 12*lambda^2*di^2 - 2*lambda^4*di^4) + 3)/ (4*lambda);
    
            W = (4*sigma^2*lambda^5)/3 * W;
    
            Wis(:,:,i) = W;
    
    
            G = zeros(3,3);
            G(1,1) = lambda^2 * di^2 + 2 * lambda + 2;
            G(1,2) = 2*(lambda*di^2+di);
            G(1,3) = di^2;
            G(2,1) = -lambda^3*di^2;
            G(2,2) = -2*(lambda^2*di^2-lambda*di-1);
            G(2,3) = 2*si - lambda*di^2;
            G(3,1) = lambda^4*di^2 - 2*lambda^3*di;
            G(3,2) = 2*(lambda^3*di^2 - 3*lambda^2*di);
            G(3,3) = lambda^2*di^2-4*lambda*di+2;
    
            G = exp(-lambda*di) / 2 * G;
    
            Gis(:,:,i) = G;
    
    
    
        end


        



    else

        R_s = 777;

    end




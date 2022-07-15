function [U_re, P_r, P_c, rhoVect, newIndices, newRawIndices] = SparseDecomp(data, learnInfo, M, L)

    %TODO: Decomp extension Matern ability.
    %TODO: Relies on KForFullR. Should be folded into existing methods.

    n = learnInfo.N;
    D = learnInfo.d;
    capN = n * D * M * L;
    tildeN = (n * (n - 1) * M * L) / 2;

    %Construct U_re.
    U_re = zeros(capN,n);
    
    for m = 1 : M
        for j = 1 : D
            for l = 1 : L
                for i1 = 1 : n
                    for i2 = 1 : n            
                        if abs(i1 - i2) > 0 && abs(i1 - i2) <= n
                            U_re((m-1)*L*n*D + (l-1)*n*D + (j-1)*n + i1, i2) = (getSpatialCoord(data, i1, j, m, l) - getSpatialCoord(data, i2, j, m, l));
                        end                    
                    end
                end
            end
        end
    end
    
    %U_re is constructed properly and ready to return.
    %Next, make the pairwise distances.   

    ds = zeros((n*(n-1))/2,1);
    entry = 1;
    for i = 1 : n
        for j = i + 1 : n
            ds(entry) = norm(getuij(data,i,j,L,M));
            entry = entry + 1;
        end
    end
    
    [ds, indices] = sort(ds);
    
    %Build P_r indexing.
    P_r = zeros(tildeN, 2);
    rowNum = 1;
    for i = 1 : n - 1
        for j = i + 1 : n
            P_r(rowNum,1) = i;
            P_r(rowNum,2) = j;
            rowNum = rowNum + 1;
        end
    end
    
    %Sort the P matrix to reflect the sorted U_s and R_s via ds.
    P_r = P_r(indices,:);
    
    %Store the needed rhos for the sparse inverse.
    rhoVect = zeros(tildeN - 1,1);
    
    for i = 1 : tildeN - 1
        rhoVect(i,1) = KForFullR(ds(i+1), ds(i));
    end

    %Now, construct the indexing array for the g4 equation.
    nextAvailable = 1;
    foldedIndexingSet = zeros(n, n-1);
    
    for i = 1 : n
    
        row = zeros(n-1,1);
    
        %i-1 edge cases + capN-i next available nums at end    
        for k = 1 : i - 1
            prevRowInt = foldedIndexingSet(k, i-1);
            row(k) = prevRowInt;
        end
    
        for k = i : n - 1
            row(k) = nextAvailable;
            nextAvailable = nextAvailable + 1;
        end
    
        foldedIndexingSet(i,:) = row;
    
    end
    
    %At this point, folded indexing set just needs to be renumbered.
    P_c = zeros(n, n-1);

    for i = 1 : n
    
        row = foldedIndexingSet(i,:);
    
        for j = 1 : n-1    
            P_c(i,j) = find(indices==row(j));    
        end
    
        P_c(i,:) = sort(P_c(i,:));
    
    end
    
    %Now P_c is constructed to specs.

    %Raw indices for permutation formatting transfer.
    pairRawIndices = zeros(n,n-1);
    
    for i = 1 : n
        pairRawIndices(i,:) = 1 : n-1;
    end
    
    %Output storage.
    newIndices = zeros(n,n-1);
    newRawIndices = zeros(n,n-1);
    
    for rowNum = 1 : n
    
        row = foldedIndexingSet(rowNum,:);
        rawRow = pairRawIndices(rowNum,:);
    
        %Get indices of each of the entries in row in the ds vector.
        ind = zeros(n-1,1);
        for i = 1 : n-1    
            target = row(i);
            for j = 1 : length(ds)    
                if isequal(indices(j), target)   
                    ind(i) = j;    
                end    
            end    
        end    
    
        %Now, sort the indices list, and reuse the sort to sort row.
        [~, thisSort] = sort(ind);
        
        %Sort row.
        row = row(thisSort);
        rawRow = rawRow(thisSort);
    
        %Save the new row.
        newIndices(rowNum,:) = row;
        newRawIndices(rowNum,:) = rawRow;
    
    
    end
    
    
    
   


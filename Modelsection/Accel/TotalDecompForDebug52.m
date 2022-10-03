function [U_s, R_s, UnsortU, ds] = TotalDecompForDebug52(data, dataA, omega, delta, n, D, M, L)


    capN = n * D * M * L;
    tildeN = (n * (n - 1) * M * L) / 2;
    
    U_s = zeros(capN,tildeN);
    
    
    %Try to build U_s with entries that are the distance pairs
    colNum = 1;
    for m = 1 : M
        for l = 1 : L
            for i = 1 : n-1
                for j = i + 1 : n
            
                    U_s(((D * i - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*i + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(dataA,i,j,l,m);
                    U_s(((D * j - 1) + (l-1)*n*D + (m-1)*n*D*L):(D*j + (l-1)*n*D + (m-1)*n*D*L),colNum) = getuij(dataA,j,i,l,m);
            
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
    
% 
%     ds
% 
% % data1 = rand(1,20)./2;      %# Sample data set 1
% % data2 = 0.3+rand(1,20)./2;  %# Sample data set 2
% hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%              'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%              'XLim',[0 1],...               %#   set the x axis limit,
%              'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
%              'Color','none');               %#   and don't use a background color
% % plot(data1,0,'r*','MarkerSize',10);  %# Plot data set 1
% % plot(data2,0,'b.','MarkerSize',10);  %# Plot data set 2
% plot(ds,0,'r*','MarkerSize',10);  %# Plot data set 1
% 
%     return
    %Swap cols
    %2 lines debug
    UnsortU = U_s;
    sortingInd = indices;
    U_s = U_s(:,indices);
    
    R_s = zeros(tildeN, tildeN);
    
    for i = 1:tildeN
        for j = 1:tildeN
    
            R_s(i,j) = KForFullRExp52(ds(i),ds(j), omega, delta,n);
    
        end
    end
    
    






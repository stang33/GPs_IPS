function rcount = obs2pdistcount(obs,edges,N)
% compute histcount for long trajectory: cut to pieces to avoid (steps*N)^2
% being too large
% size(obs)= [Nd,L] 

% (c) XXXX


[~, L] = size(obs); 

tolsize  = 1e6; pdsize = N*(N-1)/2; pieceLength = ceil(tolsize/pdsize);
if L> pieceLength
    M = ceil(L/pieceLength); 
    rcount = zeros(M,length(edges)-1);
    for m=1:M
        ind = (m-1)*pieceLength+1:min(m*pieceLength, L); 
        obstemp     = obs(:,ind);
        pdist1path  = Pairwise_distance(obstemp,N);  
        rcount(m,:) = histcounts(pdist1path,edges);
    end
    rcount = sum(rcount); 
else
    pdist1path  = Pairwise_distance(obs,N); 
    rcount      = histcounts(pdist1path,edges);
end
end

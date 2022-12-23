function K = TrueKCols(learnInfo,hyp,i,jitter)

    X = learnInfo.X;
    dN = learnInfo.d*learnInfo.N*learnInfo.order;
    L = length(X)/dN;
    d = learnInfo.d;
    N = learnInfo.N;
    K = zeros(d*N,L*d*N); % kernelmatrix
    
    
    v = learnInfo.v;
    order = learnInfo.order;
    name = learnInfo.name;
    dN1 = d*N;
    dNs = dN1;
    sub = 1:N;

    %     for j = 1:L
    %         K(1:dN1,(dN1*(j-1)+1):(dN1*j)) = k_phi_Matern_NoDeriv(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(1:2),order,name);
    %     end
    

        for j = 1:L

            %K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)) = k_phi_Matern_NoDeriv(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(1:2),order,'E',sub);
            K(1:dN1,(dN1*(j-1)+1):(dN1*j)) = k_phi_Matern_NoDeriv(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(1:2),order,'E',sub);
            
            if order == 2

                K(1:dN1,(dN1*(j-1)+1):(dN1*j)) = K(1:dN1,(dN1*(j-1)+1):(dN1*j)) + k_phi_Matern_NoDeriv(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(3:4),order,'A',sub);
               
            end
        
        end

    
    
    id = eye(dN1*L);
    K = K + jitter*id((dN1*(i-1)+1):(dN1*i),:);


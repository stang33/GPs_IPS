function K = TrueKColsALT(learnInfo,hyp,i,jitter, data, p1, p2)

    X = learnInfo.X;
    dN = learnInfo.d*learnInfo.N*learnInfo.order;
    L = length(X)/dN;
    d = learnInfo.d;
    N = learnInfo.N;
    %K = zeros(d*N,L*d*N); % kernelmatrix
    K = zeros(2,L*d*N);
    
    v = learnInfo.v;
    order = learnInfo.order;
    name = learnInfo.name;
    dN1 = d*N;
    dNs = dN1;
    sub = 1:N;

    

    ijblock = zeros(d,d);

    x1p = data(2*p1-1:2*p1, 1, 1); %TODO what about other L, M
    x2p = data(2*p2-1:2*p2, 1, 1);

    for j = 1:L

    %+1/n^2*sigma^2*exp(-omega*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*rk1x1*rk2x2';
    

        for j = 1:L

            %K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)) = k_phi_Matern_NoDeriv(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(1:2),order,'E',sub);
            K(1:2,(dN1*(j-1)+1):(dN1*j)) = k_phi_Matern_NoDeriv_ALT(X((2*(i-1)+1):(2*i)),X((2*(j-1)+1):(2*j)),learnInfo.d,v,hyp(1:2),order,'E',sub);
            
            if order == 2

                K(1:2,(dN1*(j-1)+1):(dN1*j)) = K(1:2,(dN1*(j-1)+1):(dN1*j)) + k_phi_Matern_NoDeriv_ALT(X((2*(i-1)+1):(2*i)),X((2*(j-1)+1):(2*j)),learnInfo.d,v,hyp(3:4),order,'A',sub);
               
            end
        
        end

    
    
    id = eye(dN1*L);
    K = K + jitter*id(1:2,:);


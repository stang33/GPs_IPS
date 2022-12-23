function [learnInfo] = getK_FM(learnInfo,hyp)

% input - X: training position data
%         Y: taining velocity data
%         option: 'alldata' or 'subset'

% (c) XXXX

X = learnInfo.X;
Y = learnInfo.Y;
dN = learnInfo.d*learnInfo.N*learnInfo.order;
L = length(X)/dN;
d = learnInfo.d;
N = learnInfo.N;
% CoefM = zeros(L*d*N,1); % coeficientmatrix

 Nsub = N;
 sub = 1:learnInfo.N;


K = zeros(L*d*Nsub,L*d*Nsub); % kernelmatrix

v = learnInfo.v;
order = learnInfo.order;
name = learnInfo.name;

dN1 = d*N;
dNs = d*Nsub;


if order == 2
    K_A = zeros(L*dNs,L*dNs); % kernelmatrix for phi_A
end


for i = 1:L
    [K((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)), dKs_E((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)),dKo_E((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i))] = k_phi_Matern(X((dN*(i-1)+1):(dN*i)),X((dN*(i-1)+1):(dN*i)),learnInfo.d,v,hyp(1:2),order,'E',sub);
    if order == 2
        [K_A((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)), dKs_A((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)),dKo_A((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i))] = k_phi_Matern(X((dN*(i-1)+1):(dN*i)),X((dN*(i-1)+1):(dN*i)),learnInfo.d,v,hyp(3:4),order,'A',sub);
        K((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)) = K((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i)) + K_A((dNs*(i-1)+1):(dNs*i),(dNs*(i-1)+1):(dNs*i));
    end
    for j = i+1:L
        [K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)), dKs_E((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)),dKo_E((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j))] = k_phi_Matern(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(1:2),order,'E',sub);
       
        
        if order == 2
            [K_A((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)), dKs_A((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)),dKo_A((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j))] = k_phi_Matern(X((dN*(i-1)+1):(dN*i)),X((dN*(j-1)+1):(dN*j)),learnInfo.d,v,hyp(3:4),order,'A',sub);
           
            K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)) = K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j)) + K_A((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j));
        end
        
        K((dNs*(j-1)+1):(dNs*j),(dNs*(i-1)+1):(dNs*i)) = K((dNs*(i-1)+1):(dNs*i),(dNs*(j-1)+1):(dNs*j))';

    end
end




learnInfo.K = K;


% Cholesky factorisation
[R,p] = chol(K,'lower'); % K=R'*R
if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

if strcmp('ODS',name)
    Ym = zeros(dNs*L,1);
    Ydk = zeros(dNs*L,1);
    YdPa = zeros(dNs*L,1);
    logk = hyp(4); % parameter in control term

    YdPb = 0*zeros(dNs*L,1);
    YdPc = 0*zeros(dNs*L,1);
    Pa = hyp(5); % parameter in control term
    Pb = hyp(6); % parameter in control term
    Pc = hyp(7); % parameter in control term
%     for i = 1:L
%         Ym((dN1*(i-1)+1):(dN1*(i-1)+d)) = Y((dN1*(i-1)+1):(dN1*(i-1)+d)) + exp(logk)*(X((dN1*(i-1)+1):(dN1*(i-1)+d)) - Pa);
%         Ym((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) = Y((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) + exp(logk)*(X((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) - Pb);
%         Ym((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) = Y((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) + exp(logk)*(X((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) - Pc);
%         Ydk((dN1*(i-1)+1):(dN1*(i-1)+d)) = (X((dN1*(i-1)+1):(dN1*(i-1)+d)) - Pa)*exp(logk);
%         Ydk((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) = (X((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) - Pb)*exp(logk);
%         Ydk((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) = (X((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) - Pc)*exp(logk);
%         YdPa((dN1*(i-1)+1):(dN1*(i-1)+d)) = -exp(logk);
%         YdPb((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) = -exp(logk);
%         YdPc((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) = -exp(logk);
%     end

    for i = 1:L
        Ym((dNs*(i-1)+1):(dNs*(i-1)+d)) = Y((dNs*(i-1)+1):(dNs*(i-1)+d)) + exp(logk)*(X((dNs*(i-1)+1):(dNs*(i-1)+d)) - Pa);
        Ym((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) = Y((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) + exp(logk)*(X((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) - Pb);
        Ym((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) = Y((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) + exp(logk)*(X((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) - Pc);
        Ydk((dNs*(i-1)+1):(dNs*(i-1)+d)) = (X((dNs*(i-1)+1):(dNs*(i-1)+d)) - Pa)*exp(logk);
        Ydk((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) = (X((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) - Pb)*exp(logk);
        Ydk((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) = (X((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) - Pc)*exp(logk);
        YdPa((dNs*(i-1)+1):(dNs*(i-1)+d)) = -exp(logk);
        YdPb((dNs*(i-1)+d+1):(dNs*(i-1)+2*d)) = -exp(logk);
        YdPc((dNs*(i-1)+2*d+1):(dNs*(i-1)+3*d)) = -exp(logk);
    end

    learnInfo.Ym = Ym;
    alpha = R'\(R\Ym);
    fval = 0.5*Ym'*alpha + sum(log(diag(R))) + log(2*pi)*(dN1*L)/2;

    dfval = 0*hyp;

    Q =  R'\(R\eye(L*dN1)) - alpha*alpha';
    dfval(1) = sum(sum(Q.*dKs_E))/2;
    dfval(2) = sum(sum(Q.*dKo_E))/2;

    if ~isnan(hyp(3))
        dfval(3) = 2*exp(lognoise)^2*trace(Q)/2;
    else
        dfval(3) = 0;
    end

    dfval(4) = alpha'* Ydk;
    dfval(5) = alpha'* YdPa;
    dfval(6) = alpha'* YdPb;
    dfval(7) = alpha'* YdPc;


elseif strcmp('CSF',name)
    Ym = zeros(dNs*L,1);
    Yda = zeros(dNs*L,1);
    Ydb = zeros(dNs*L,1);
    loga = hyp(6); % parameter in noncollective force
    logb = hyp(7);
    for i = 1:L
        for j = 1:Nsub
            Ym((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j))) - (exp(loga) * X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*(1- norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb)));
            Yda((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = - X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*(1- norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb))*exp(loga);
            Ydb((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = exp(loga) * X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb)*log(norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2))*exp(logb);
        end
    end
    learnInfo.Ym = Ym;
    alpha = R'\(R\Ym);
    fval = 0.5*Ym'*alpha + sum(log(diag(R))) + log(2*pi)*(dN1*L)/2;

    dfval = 0*hyp;

    Q =  R'\(R\eye(L*dNs)) - alpha*alpha';
    dfval(1) = sum(sum(Q.*dKs_E))/2;
    dfval(2) = sum(sum(Q.*dKo_E))/2;
    dfval(3) = sum(sum(Q.*dKs_A))/2;
    dfval(4) = sum(sum(Q.*dKo_A))/2;

    if ~isnan(hyp(5))
        dfval(5) = 2*exp(lognoise)^2*trace(Q)/2;
    else
        dfval(5) = 0;
    end

    dfval(6) = alpha'* Yda;
    dfval(7) = alpha'* Ydb;


elseif strcmp('FM',name)



    Ym = zeros(dNs*L,1);
    Yda = zeros(dNs*L,1);
    Ydb = zeros(dNs*L,1);
    loga = hyp(6); % parameter in noncollective force
    logb = hyp(7);

    for i = 1:L
        for j = 1:Nsub
            Ym((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j))) - (exp(loga) - exp(logb)*sum(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))).^2))*X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)));
            Yda((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = - X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*exp(loga);
            Ydb((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = sum(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))).^2)*X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*exp(logb);

        end
    end
    learnInfo.Ym = Ym;
  
  

    
else
    Ym = zeros(dNs*L,1);
    for i = 1:L
        for j = 1:Nsub
            Ym((dNs*(i-1)+1+d*(j-1)):(dNs*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j)));
        end
    end

    learnInfo.Ym = Ym;
    alpha = R'\(R\Ym);
    fval = 0.5*Ym'*alpha + sum(log(diag(R))) + log(2*pi)*(dNs*L)/2;

    dfval = 0*hyp;

    Q =  R'\(R\eye(L*dNs)) - alpha*alpha';
    dfval(1) = sum(sum(Q.*dKs_E))/2;
    dfval(2) = sum(sum(Q.*dKo_E))/2;
    dfval(3) = sum(sum(Q.*dKs_A))/2;
    dfval(4) = sum(sum(Q.*dKo_A))/2;

    if ~isnan(hyp(5))
        dfval(5) = 2*exp(lognoise)^2*trace(Q)/2;
    else
        dfval(5) = 0;
    end

end

% learnInfo.CoefM = pinv(K)*Ym;

end




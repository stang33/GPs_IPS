function [mf] = predict_ode(xstar,learnInfo)
% function [mf, covf,covfxi] = predict_ode(xstar,learnInfo)

% Input xstat: dN x 1 vector 
% (c) XXXX

X = learnInfo.X;
Y = learnInfo.Y;
N = learnInfo.N;
d = learnInfo.d;
order = learnInfo.order;
name = learnInfo.name;
dN = d*N*order;
L = learnInfo.L;
K = learnInfo.K;
hyp = learnInfo.hyp;
CoefM = learnInfo.CoefM;


if strcmp('subset',learnInfo.option)
    Nsub = learnInfo.Nsub;
    sub = learnInfo.sub;
else
    Nsub = N;
    sub = 1:learnInfo.N;
end

dNs = d*Nsub;
dN1= d*N;
K_s_E = zeros(dN1,L*dNs); %covariance matrix between new point (xstar) and data (X)


for j = 1:L
    [K_s_E(:,(dNs*(j-1)+1):(dNs*j)),~,~] = k_phi_Matern(xstar,X((dN*(j-1)+1):(dN*j)),learnInfo.d,learnInfo.v,hyp(1:2),order,'E',1:N,sub);
end

% [K_ss_E,~,~] = k_phi_Matern(xstar,xstar,learnInfo.d,learnInfo.v,hyp(1:2),order,'E');

if order == 2
    K_s_A = zeros(dN1,L*dNs);
    for j = 1:L
        [K_s_A(:,(dNs*(j-1)+1):(dNs*j)),~,~] = k_phi_Matern(xstar,X((dN*(j-1)+1):(dN*j)),learnInfo.d,learnInfo.v,hyp(3:4),order,'A',1:N,sub);
    end
    % [K_ss_A,~,~] = k_phi_Matern(xstar,xstar,learnInfo.d,learnInfo.v,hyp(3:4),order,'A');
end


if order == 1
    if strcmp(name,'ODS')
        % Ym = learnInfo.Ym;
        % mf = K_s_E*pinv(K)*Ym;
        mf = K_s_E*CoefM;
        logk = hyp(4);
        Pa = hyp(5);
        mf(1:d,1) = mf(1:d,1) - exp(logk)*(xstar(1:d) - Pa);
        Pb = hyp(6);
        Pc = hyp(7);
        mf(d+1:2*d,1) = mf(d+1:2*d,1) - exp(logk)*(xstar(d+1:2*d) - Pb);
        mf(2*d+1:3*d,1) = mf(2*d+1:3*d,1) - exp(logk)*(xstar(2*d+1:3*d) - Pc);
         
        % covf = diag(K_ss - K_s_E*pinv(K)*K_s_E');
    else    
    % mf = K_s_E*pinv(K)*Y;
    mf = K_s_E*CoefM;
    % covf = diag(K_ss_E - K_s_E*pinv(K)*K_s_E');
    end
end

if order == 2 
    if strcmp(name,'CSF')
        % Ym = learnInfo.Ym;
        mf(1:dN/order,1)=xstar(dN/order+1:end);
        % mfphi = K_s_E*pinv(K)*Ym + K_s_A*pinv(K)*Ym;
        mfphi = K_s_E*CoefM + K_s_A*CoefM;
        loga = hyp(6);
        logb = hyp(7);
        for j = 1:N
            mf(dN/order+1+d*(j-1):dN/order+d*j,1) = exp(loga)*xstar(dN/order+1+d*(j-1):dN/order+d*j)*(1-norm(xstar(dN/order+1+d*(j-1):dN/order+d*j),2)^exp(logb)) + mfphi(1+d*(j-1):d*j,1);
        end
        % covf = diag(K_ss_E - K_s_E*pinv(K)*K_s_E') + diag(K_ss_A - K_s_A*pinv(K)*K_s_A');
        
	elseif strcmp(name,'FM')
        % Ym = learnInfo.Ym;
        mf(1:dN/order,1)=xstar(dN/order+1:end);
        % mfphi = K_s_E*pinv(K)*Ym + K_s_A*pinv(K)*Ym;
        mfphi = K_s_E*CoefM + K_s_A*CoefM;
        loga = hyp(6);
        logb = hyp(7);
        for j = 1:N
            mf(dN/order+1+d*(j-1):dN/order+d*j,1) = (exp(loga) - exp(logb)*sum(xstar(dN/order+1+d*(j-1):dN/order+d*j).^2))*xstar(dN/order+1+d*(j-1):dN/order+d*j) + mfphi(1+d*(j-1):d*j,1);
        end
        % covf = diag(K_ss_E - K_s_E*pinv(K)*K_s_E') + diag(K_ss_A - K_s_A*pinv(K)*K_s_A');
        
        
    else
        % Ym = learnInfo.Ym;
        mf(1:dN/order,1)=xstar(dN/order+1:end);
        % mfphi = K_s_E*pinv(K)*Ym + K_s_A*pinv(K)*Ym;
        mfphi = K_s_E*CoefM + K_s_A*CoefM;
        for j = 1:N
            mf(dN/order+1+d*(j-1):dN/order+d*j,1) = mfphi(1+d*(j-1):d*j,1);
        end
        % covf = diag(K_ss_E - K_s_E*pinv(K)*K_s_E') + diag(K_ss_A - K_s_A*pinv(K)*K_s_A');
        
    end
end


end
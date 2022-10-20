function [Kr] = K_r(r,X,learnInfo, kernel_type)
% Input: r:  1 x 1 
%        X:  dNL x 1  position data
%        kernel_type: kernel for 'E' or 'A'
% Output: K_r: 1 x dN covariance function

% (c) XXXX

N = learnInfo.N;
d = learnInfo.d;
hyp = learnInfo.hyp;
%Kr = zeros(length(r),length(X)/learnInfo.order);

if strcmp('subset',learnInfo.option)
    Nsub = learnInfo.Nsub;
    sub = learnInfo.sub;
else
    Nsub = N;
    sub = 1:learnInfo.N;
end



v = learnInfo.v;
order = learnInfo.order;
d = learnInfo.d;
N = learnInfo.N;
dN= d*N*order;


if order == 1 || kernel_type == 'E'
    sigma = exp(hyp(1));
    omega = exp(hyp(2));
else
    sigma = exp(hyp(3));
    omega = exp(hyp(4));
end


for j=1:length(r)
    for l=1:length(X)/dN
        X_state = reshape(X((l-1)*dN+1:l*dN,:),d,[]); % d x (N xorder)
        for i =1: Nsub
            r_temp = X_state(:,1:N)-repmat(X_state(:,sub(i)),1,N); % [x_1-x_i; x_2-x_i; \cdots, x_N-x_i]  d x N
            r_temp_norm = sqrt(sum(r_temp.^2,1)); % 1 x N
            if order ==2 && kernel_type == 'A'
            r_temp = X_state(:,N+1:end)-repmat(X_state(:,sub(i)+N),1,N); % [v_1-v_i; v_2-v_i; \cdots, v_N-v_i]  d x N
            end

            Kr(j,(l-1)*d*Nsub+1+(i-1)*d:(l-1)*d*Nsub+i*d) = (r_temp*cov_Matern(r_temp_norm,r(j),sigma,omega,v)'/N)';                         
        end
    end
end


end

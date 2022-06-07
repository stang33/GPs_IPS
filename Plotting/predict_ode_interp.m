function [mf] = predict_ode_interp(xstar,learnInfo)

% Input xstat: dN x 1 vector 
% (c) XXXX

X = learnInfo.X;
Y = learnInfo.Y;
N = learnInfo.N;
d = learnInfo.d;
order = learnInfo.order;
name = learnInfo.name;
dN = d*N*order;
hyp = learnInfo.hyp;



x1p = reshape(xstar(1:N*d),d,[]); %d-by-n matrix
if order == 2
    x1v =reshape(xstar(N*d+1:2*N*d),d,[]); %d-by-n matrix
end


if order ==1
    aa= zeros(d,N);
    for i=1:N
        temp      = x1p- repmat(x1p(:,i),1,N);  % d x N [x_1-x_i,\cdots, x_N-x_i]
        DD        = sqrt( sum(temp.^2,1) ); % 1 x N
        DD(DD==0) = 1;  %  this step is to avoid 0 x inf= NAN
        aa(:,i)   = interp1(learnInfo.r,learnInfo.phi',DD)'/N;   % d x 1
    end
    mf = reshape(aa,[],1); % dN x 1   
    if strcmp(name,'ODS')
        logk = hyp(4);
        Pa = hyp(5);
        mf(1:d,1) = mf(1:d,1) - exp(logk)*(xstar(1:d) - Pa);
        Pb = hyp(6);
        Pc = hyp(7);
        mf(d+1:2*d,1) = mf(d+1:2*d,1) - exp(logk)*(xstar(d+1:2*d) - Pb);
        mf(2*d+1:3*d,1) = mf(2*d+1:3*d,1) - exp(logk)*(xstar(2*d+1:3*d) - Pc);
    end
end

if order ==2 
    aa= zeros(d,N);

    if strcmp(name,'CSF')
        for i=1:N
            temp_position      = x1p(:,[1:i-1 i+1:N])- repmat(x1p(:,i),1,N-1);  % d x N-1 [x_1-x_i,\cdots,x_{i-1}-x_i,x_{i+1}-x_i,\cdot,x_N-x_i]
            DD        = sqrt( sum(temp_position.^2,1) ); % 1 x N-1
            temp_velocity      = x1v(:,[1:i-1 i+1:N])- repmat(x1v(:,i),1,N-1);  % d x N-1 [v_1-v_i,\cdots,v_{i-1}-v_i,v_{i+1}-v_i,\cdot, v_N-v_i]
            aa(:,i)   = temp_velocity * interp1(learnInfo.r,learnInfo.phi',DD)'/N;   % 1 x d
        end
        mf(1:dN/order,1)=xstar(dN/order+1:end);
        mfphi = reshape(aa,[],1); % dN x 1  
        loga = hyp(4);
        logb = hyp(5);
        for j = 1:N
            mf(dN/order+1+d*(j-1):dN/order+d*j,1) = exp(loga)*xstar(dN/order+1+d*(j-1):dN/order+d*j)*(1-norm(xstar(dN/order+1+d*(j-1):dN/order+d*j),2)^exp(logb)) + mfphi(1+d*(j-1):d*j,1);
        end
        % covf = diag(K_ss - K_s*invK*K_s');
        
    elseif strcmp(name,'FM')
        for i=1:N
            temp_position      = x1p(:,[1:i-1 i+1:N])- repmat(x1p(:,i),1,N-1);  % d x N-1 [x_1-x_i,\cdots,x_{i-1}-x_i,x_{i+1}-x_i,\cdot,x_N-x_i]
            DD        = sqrt( sum(temp_position.^2,1) ); % 1 x N-1
            aa(:,i)   = temp_position * interp1(learnInfo.r,learnInfo.phi',DD)'/N;   % 1 x d
        end
        mf(1:dN/order,1)=xstar(dN/order+1:end);
        mfphi = reshape(aa,[],1); % dN x 1  
        
        loga = hyp(4);
        logb = hyp(5);
        for j = 1:N
            mf(dN/order+1+d*(j-1):dN/order+d*j,1) = (exp(loga) - exp(logb)*sum(xstar(dN/order+1+d*(j-1):dN/order+d*j).^2))*xstar(dN/order+1+d*(j-1):dN/order+d*j) + mfphi(1+d*(j-1):d*j,1);
        end
        %covf = diag(K_ss - K_s*invK*K_s');
        
    end
end


end
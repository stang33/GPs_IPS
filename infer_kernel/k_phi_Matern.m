function [k, dks, dko] = k_phi_Matern(x1,x2,d,v,hyp,order,kernel_type,sub1,sub2)
%X,Y are dN dimension vectors contain values of f_phi
% input  -   x1: dN dimension vector
%            x2: dN dimension vector 
%             d: dimension of position vector for a agent
%             v: 1/2,3/2,5/2,7/2
%           hyp: hyperparameter for the convariance function
%         order: order of the dynamic system
%   kernel_type: kernel for 'E' or 'A'
%           sub: the subset of data we use
%
% k(x,z) = max(1-r,0)^(j+v) * f(r,j) with j = floor(D/2)+v+1
%
% where r is the Mahalanobis distance sqrt(maha(x,z)). The hyperparameters are:
%
% hyp = [ hyp_maha ]

% (c) XXXX


%  For second order system:     X_train consists of [x_1,\cdots,x_N,v_1,\cdots,v_N]
%                               Y_train consistss of [\dot v_1,..., \dot v_N]
%
%
%  for second order system with xi X_train consists of [x_1,\cdots,x_N,v_1,\cdots,v_N,xi_1,...,xi_N]

logsigma = hyp(1);
logomega = hyp(2);
sigma = exp(hyp(1));
omega = exp(hyp(2));

n = size(x1,1)/(d*order);

if nargin<8, sub1 = 1:n; sub2 = 1:n; end
if nargin == 8, sub2 = sub1; end


if order == 2
    x1v =reshape(x1(n*d+1:2*n*d),d,[]); %d-by-n matrix
    x2v =reshape(x2(n*d+1:2*n*d),d,[]); %d-by-n matrix
    rx1v = zeros(d*n,n);   %dn-by-n matrix with directed distance of X1V
    rx2v = zeros(d*n,n);   %dn-by-n matrix with directed distance of X2V
end

x1p = reshape(x1(1:n*d),d,[]); %d-by-n matrix
x2p = reshape(x2(1:n*d),d,[]); %d-by-n matrix
rx1p = zeros(d*n,n);   %dn-by-n matrix with directed distance of X1P
rx2p = zeros(d*n,n);   %dn-by-n matrix with directed distance of X2P


for i = 1:n
    for j = 1:i-1
        
        rx1p((d*(i-1)+1):d*i,j) = x1p(:,i)-x1p(:,j);  
        rx1p((d*(j-1)+1):d*j,i) = -rx1p((d*(i-1)+1):d*i,j);
        rx2p((d*(i-1)+1):d*i,j) = x2p(:,i)-x2p(:,j);  
        rx2p((d*(j-1)+1):d*j,i) = -rx2p((d*(i-1)+1):d*i,j);
        
        if order ==2 && kernel_type == 'A' % e.g. for cuker smale and phototaxi
        rx1v((d*(i-1)+1):d*i,j) = x1v(:,i)-x1v(:,j);  
        rx1v((d*(j-1)+1):d*j,i) = -rx1v((d*(i-1)+1):d*i,j);
        rx2v((d*(i-1)+1):d*i,j) = x2v(:,i)-x2v(:,j);  
        rx2v((d*(j-1)+1):d*j,i) = -rx2v((d*(i-1)+1):d*i,j);
        end
    end
end


norm_rx1 = zeros(n,n);   %n-by-n matrix with distance of X 
norm_rx2 = zeros(n,n);   %n-by-n matrix with distance of Y 
for i = 1:n
    for j = 1:i-1
        norm_rx1(i,j) = norm(rx1p((d*(i-1)+1):d*i,j));
        norm_rx1(j,i) = norm_rx1(i,j);
        norm_rx2(i,j) = norm(rx2p((d*(i-1)+1):d*i,j));
        norm_rx2(j,i) = norm_rx2(i,j);
    end
end


Nsub1 = length(sub1);
Nsub2 = length(sub2);
k = zeros(Nsub1*d,Nsub2*d); %kernel matrix
dks = zeros(Nsub1*d,Nsub2*d); %derivative of kernel matrix w.r.t. logsigma
dko = zeros(Nsub1*d,Nsub2*d); %derivative of kernel matrix w.r.t. logomega

switch v
    case 1/2
        for i = 1:Nsub1
            for j = 1:Nsub2
                for k1 = [1:sub1(i)-1,sub1(i)+1:n]
                    for k2 = [1:sub2(j)-1,sub2(j)+1:n]
                        
                        if kernel_type == 'A'
                            rk1x1 = rx1v((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2v((d*(k2-1)+1):d*k2,sub2(j));
                        else
                            rk1x1 = rx1p((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2p((d*(k2-1)+1):d*k2,sub2(j));
                        end
                            k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+1/n^2*sigma^2*exp(-omega*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*rk1x1*rk2x2';
                            dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+1/n^2*sigma^2*exp(-omega*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*rk1x1*rk2x2'*2;
                            dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+1/n^2*sigma^2*exp(-omega*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*rk1x1*rk2x2'*(-omega)*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j)));

                    end
                end
            end
        end
        
    case 3/2
        
        for i = 1:Nsub1
            for j = 1:Nsub2
                for k1 = [1:sub1(i)-1,sub1(i)+1:n]
                    for k2 = [1:sub2(j)-1,sub2(j)+1:n]
                         
                        if kernel_type == 'A'
                            rk1x1 = rx1v((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2v((d*(k2-1)+1):d*k2,sub2(j));
                        else
                            rk1x1 = rx1p((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2p((d*(k2-1)+1):d*k2,sub2(j));
                        end
                            
                        f_val_temp = 1/n^2*exp(2*logsigma-exp(logomega)*sqrt(3)*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*(1+sqrt(3)*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j)))*exp(logomega))*rk1x1*rk2x2';
                        k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+f_val_temp;
                        dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+2*f_val_temp;
                        o_temp_vec = 1/n^2*exp(2*logsigma)*exp(-exp(logomega)*sqrt(3)*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j))))*rk1x1*rk2x2';
                        o_temp_factor = sqrt(3)* abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j)))*exp(logomega);
                        dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))-o_temp_vec*o_temp_factor^2;
                    end
                end
            end
        end
        
      case 5/2
          sigma=exp(logsigma);
          omega=exp(logomega);
        
        for i = 1:Nsub1
            for j = 1:Nsub2
                for k1 = [1:sub1(i)-1,sub1(i)+1:n]
                    for k2 = [1:sub2(j)-1,sub2(j)+1:n]
                        
                        if kernel_type == 'A'
                            rk1x1 = rx1v((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2v((d*(k2-1)+1):d*k2,sub2(j));
                        else
                            rk1x1 = rx1p((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2p((d*(k2-1)+1):d*k2,sub2(j));
                        end
                            
                        r_temp=abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j)))*sqrt(5)*omega;
                        f_val_temp = 1/n^2*sigma^2*exp(-r_temp)*(1+r_temp+r_temp^2/3)*rk1x1*rk2x2';
                        k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+f_val_temp;
                        dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+2*f_val_temp;
                        o_temp_vec = 1/n^2*sigma^2*exp(-r_temp)*rk1x1*rk2x2';
                        dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))-o_temp_vec*r_temp^2*(1+r_temp)/3;
                    end
                end
            end
        end
        
        case 7/2
        omega =exp(logomega);
        sigma=exp(logsigma);
        for i = 1:Nsub1
            for j = 1:Nsub2
                for k1 = [1:sub1(i)-1,sub1(i)+1:n]
                    for k2 = [1:sub2(j)-1,sub2(j)+1:n]
                        
                        if kernel_type == 'A'
                            rk1x1 = rx1v((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2v((d*(k2-1)+1):d*k2,sub2(j));
                        else
                            rk1x1 = rx1p((d*(k1-1)+1):d*k1,sub1(i));
                            rk2x2 = rx2p((d*(k2-1)+1):d*k2,sub2(j));
                        end
                            
                        r_temp=sqrt(7)*abs(norm_rx1(k1,sub1(i))-norm_rx2(k2,sub2(j)))*omega;
                        f_val_temp = 1/n^2*sigma^2*exp(-r_temp)*(1+r_temp+0.4*r_temp^2+r_temp^3/15)*rk1x1*rk2x2';
                        k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = k((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+f_val_temp;
                        dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dks((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))+2*f_val_temp;
                        o_temp_vec = 1/n^2*sigma^2*exp(-r_temp)*rk1x1*rk2x2';
                        dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j)) = dko((d*(i-1)+1):(d*i),(d*(j-1)+1):(d*j))-o_temp_vec*r_temp^2*(1+r_temp+r_temp^2/3)/5;
                    end
                end
            end
        end
      
       
        
end

    
    








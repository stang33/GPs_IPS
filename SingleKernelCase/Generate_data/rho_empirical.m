function rho = rho_empirical(obs_data,sysInfo,obsInfo,saveON,plotON)
% Estimate empirical rho, density of pairwise distance, from data 
% Input :  
%        the trajectory data obs_data: Nd x L x  M
%       
% Output: 
%        I: Rmin:h:Rmax
%        density: emprical rho on I
%        p1: pairewise distance of vq1

% (c) XXXX

[~,~,M] = size(obs_data); % L = # of time steps; M = number of trajectories
N        = sysInfo.N;        % number of agents
d        = sysInfo.d;

edges  = obsInfo.rho_T_histedges;  
rcount = zeros(M,length(edges)-1);     % histcount of rho for each path 

binsize = edges(2:end)-edges(1:end-1); 

parfor i = 1:M  % can use parfor later
    rcount(i,:) = obs2pdistcount(obs_data(1:N*d,:,i),edges,N);
  
end


if M>1
    temp   = sum(rcount); 
else
    temp   = rcount;
end

rdistr = temp/sum(temp);               % prob r in a bin  
rdens  = rdistr./binsize;              % density rho  = prob/ size(bin)

imin      = find(rdens>0,1,'first');
imax      = find(rdens>0,1,'last');
rhoSupp   = rdens(imin:imax);   % rho on support
edgesSupp = edges(imin:imax);   % edges on support of rho

rho.edges     = edges;      
rho.rdens     = rdens; 
rho.edgesSupp = edgesSupp;  
rho.rhoSupp   = rhoSupp;  

if nargin<5 || saveON ==1 
    if ~exist('outputs','dir'), mkdir('outputs'); end
    filename_temp=strcat('outputs/',sysInfo.name,'rho_empirical.mat');
    save(filename_temp, 'rho', 'sysInfo', 'obsInfo');
end

if nargin>3 && plotON ==1
    titl = strcat('The empirical rho with M = ',num2str(M));
    figure; bar(rho.edgesSupp,rho.rhoSupp); title(titl); 
    filename_temp=strcat('outputs/fig_',sysInfo.name,'_rho_empirical.eps');
    print(filename_temp,'-depsc');
end


end

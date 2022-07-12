function [e1, e2] = errornorms_phis(sysInfo,obsInfo,learnInfo,range)

% compute L-infinity and L2-rhoT norms for phis

% (c) M. Zhong, S.Tang, Jinchao(JHU)

%% load the parameters
rho_emp = learnInfo.rhoLT;
one_knot                =  rho_emp.edges;% find out the knot vector
if nargin<4,  range   = [one_knot(1), one_knot(end)];  end                               % plot the uniform learning result first; use the knot vectors as the range
r           = one_knot(1):(one_knot(2)-one_knot(1)):one_knot(end);                                % refine this knot vector, so that each sub-interval generated a knot has at least 7 interior points, so 3 levels of refinement
%r(r<range(1) | r>range(2))  = [];
%r                           = r(r>=0 & r<=rho_emp.edgesSupp(end) );

%% for phi in the state variable 
phi                       = sysInfo.phi{1}(r);
[phi_mean,~]=phi_pos(r,learnInfo);


edges               = obsInfo.rho_T_histedges;  % Estimated \rho's
edges_idxs_fine          = find(one_knot(1) <= edges & edges<one_knot(end));
%centers_fine             = (edges(edges_idxs_fine(1):edges_idxs_fine(end)-1) + edges(edges_idxs_fine(1)+1:edges_idxs_fine(end)))/2;
% downsampling of hist data by 50
edges_idxs =edges_idxs_fine(1:5:end);
%centers =centers_fine(1:20:end);
centers =edges(edges_idxs);
%histdata1            = rho_emp.rdens(edges_idxs(1:end));                    % this is the "true" \rhoLT from many MC simulations

%% compute the error norms for phis
Linfinity_idxs          = find(range(1) <= edges & edges<range(2));
e1 = max(abs(phi(Linfinity_idxs)-phi_mean(Linfinity_idxs)'))./max(abs(phi(Linfinity_idxs)));     %L-infinity norm
e2 = sqrt(sum((phi(2:end)-phi_mean(2:end)').^2.*rho_emp.rdens*(one_knot(2)-one_knot(1))))./sqrt(sum((phi(2:end)-0*phi_mean(2:end)').^2.*rho_emp.rdens*(one_knot(2)-one_knot(1))));  %L2rhoT norm
%e2 = sqrt(sum((phi(edges_idxs)-phi_mean(edges_idxs)').^2.*rho_emp.rdens(edges_idxs)*(centers(2)-centers(1)))); % approximated L2rhoT norm



end

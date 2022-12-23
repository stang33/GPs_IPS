function [ ind,dists ] = FPS( X, k )
% Farthest Point Sampling(FPS)
% sample k points from X using FPS
% X is n-by-d, n=#pts
% Difeng Cai 06-15-2021

[n,d] = size(X);
ind = ones(1,k);
S = zeros(k,d);
dists = zeros(1,k);

% initialization
X1 = X; JJ = 1:n;

% ind(1) = 1;
% S(1,:) = X1(ind(1),:);

% use initial point closest to mean(X)
[~,ind(1)] = pdist2( X,mean(X), 'Euclidean', 'Smallest', 1);
S(1,:) = X1(ind(1),:);

JJ = setdiff(JJ, ind(1));
X1 = X(JJ,:);

for i = 2:k
    % for each x in X1, find dist(x,S)
    dd = pdist2( S(1:i-1,:),X1, 'Euclidean', 'Smallest', 1);

    % find the largest distance in dd, index in JJ (or X1)
    [thedist,foo] = max(dd);

    % index of the point in X
    ind(i) = JJ(foo);
    S(i,:) = X1(foo,:);
    dists(i) = thedist;

    JJ = setdiff(JJ, ind(i));
    X1 = X(JJ,:);
end
function w = stlsq_solver(Phi, y, lambda, maxit)
% Sequential thresholded least squares
%
% Inputs:
%   Phi    : design matrix
%   y      : target
%   lambda : threshold
%   maxit  : number of iterations
%
% Output:
%   w      : sparse coefficient vector

    if nargin < 4
        maxit = 10;
    end

    % initial least-squares
    w = Phi \ y;

    for it = 1:maxit
        small = abs(w) < lambda;
        if all(~small)
            break;
        end
        w(small) = 0;

        big = ~small;
        if ~any(big)
            break;
        end

        w(big) = Phi(:,big) \ y;
    end
end
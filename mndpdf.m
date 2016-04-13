function p = mndpdf(x, X)

% mndpdf
%
% Probability density function of a multivariate normal distribution.
% Returns the likelihood corresponding to each vector of |x|.
%
% Inputs:
%
% x   Matrix whose columns are the vector locations of each point at which
%     to evalaute the PDF (nx-by-n)
% X   Covariance matrix of the distribution (strictly positive-definite)
%     (nx-by-nx)
%
% Outputs:
%
% p   Probability density at each point (1-by-n)
% 

% Copyright 2016 An Uncommon Lab

%#codegen

    % Dimensions
    [nd, n] = size(x);
    
    % Scale factor
    invR = inv(X);
    s    = 1/sqrt((2*pi)^nd * det(X));
    
    % PDF at each point.
    p = zeros(1, n);
    for k = 1:n
        p(k) = s * exp(-0.5 * x(:,k).' * invR * x(:,k)); %#ok<MINV>
    end
    
end % mndpdf

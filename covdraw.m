function c = covdraw(C, n)

% covdraw
%
% Draw n random samples from a zero-mean Gaussian distribution with 
% covariance C.
% 
% Inputs:
%
% C    Covariance matrix (nC-by-nC)
% n    Number of samples to draw (default is 1)
%
% Outputs:
%
% c    Draws from C (nC-by-n)
% 

% Copyright 2016 An Uncommon Lab

    % Defaults
    if nargin < 2, n = 1; end;

    % Use singular value decomposition instead of Cholesky, because
    % Cholesky chokes when the matrix is positive semi-definite, and we
    % want to be able to handle positive semi-definite matrices.
    [u, s]   = svd(C);
    s(s < 0) = 0; % Nothing should be negative, but just in case...
    s        = sqrt(s);
    
    % At this point, C == (u * s) * (u * s).', so (u * s) is the Cholesky
    % factor (if C was positive definite).
    
    % Construct the draws.
    c = u * s * randn(size(C, 1), n);

end % covdraw

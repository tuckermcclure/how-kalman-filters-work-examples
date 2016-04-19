function ell = ellipse(P, x, n, indices)
    
% ellipse
% 
% Create the points of the 3-sigma ellipse indicated by covariance matrix P
% centered at x. The ellipse will be 2D and will correspond to P(1:2, 1:2)
% or P(indices, indices) when indices are provided.
%
% An optional number of points, n, can be provided.
% 
% Inputs:
%
% P        Covariance matrix
% x        Center for ellipse
% n        Number of points to use in the ellipse (default is 100)
% indices  Indices to use for 2D ellise (default is 1:2)
%
% Outputs:
%
% ell  Points of the ellipse (2-by-n)
% 

% Copyright 2016 An Uncommon Lab

    % Set default number of points and indices.
    if nargin < 3 || isempty(n),       n = 100;       end;
    if nargin < 4 || isempty(indices), indices = 1:2; end;
    
    % Create a unit circle.
    theta = linspace(0, 2*pi, n);
    ell   = [cos(theta); sin(theta)];
    
    % Scale the unit circle appropriately and add on x.
    C   = chol(P(indices,indices), 'lower');
    ell = 3 * C * ell;                       % ellipse about (0, 0)
    ell = bsxfun(@plus, ell, x(indices));    % ellipse about x

end % ellipse

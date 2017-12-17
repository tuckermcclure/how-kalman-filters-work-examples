function [C, v] = randcov(n, eig_min, eig_max, type)

% randcov
%
% Create a random |n|-by-|n| covariance matrix with eigenvalues between the 
% given bounds.
% 
% Inputs:
% 
% n         Dimension of the covariance
% eig_min   Minimum allowable eigenvalue, a scalar or 1-by-|n| vector of
%           minima; should be >= 0 for a covariance matrix (default 0)
% eig_max   Maximum allowable eigenvalue, a scalar or 1-by-|n| vector of
%           maxima (default 1)
% type      |'cov'| for a covariance matrix, |'sqrt'| for a matrix square
%           root of the covariance
% 
% Outputs:
% 
% C         |n|-by-|n| covariance matrix with random eigenvalues
% v         1-by-|n| eigenvalues
% 
% This function comes from *kf, a design tool for state estimators:
% 
%    http://www.anuncommonlab.com/starkf/
% 
% Example:
%
% Create a covariance matrix, make random draws, and see that their
% covariance is approximate as requested.
% 
% % Generate the covariance matrix
% C = randcov(5)
% 
% % Draw a bunch of random samples from a binomial distribution.
% d = covdraw(C, 100000);
% 
% % Find the covariance of the draws; note that it looks like C.
% C_empirical = cov(d.')
% 
% [Command Window outputs]
% | C =
% |     0.2547    0.0298    0.1850   -0.1950   -0.0339
% |     0.0298    0.2380   -0.0108    0.0758    0.1653
% |     0.1850   -0.0108    0.3915   -0.3330   -0.0362
% |    -0.1950    0.0758   -0.3330    0.4366    0.0326
% |    -0.0339    0.1653   -0.0362    0.0326    0.2857
% | C_empirical =
% |     0.2559    0.0296    0.1852   -0.1950   -0.0346
% |     0.0296    0.2391   -0.0108    0.0755    0.1651
% |     0.1852   -0.0108    0.3917   -0.3324   -0.0360
% |    -0.1950    0.0755   -0.3324    0.4352    0.0319
% |    -0.0346    0.1651   -0.0360    0.0319    0.2858
% 
% Create a "matrix square root" of a covariance matrix with eigenvalues 
% between 1 and 4 and get the generated eigenvalues too.
% 
% [sqC, v] = randcov(4, 1, 4, 'sqrt');
%  
% % Calculate the full covariance matrix and get its eigenvalues.
% C = sqC * sqC.';
% eig(C).'
%  
% % Compare to the eigenvalues used to generate sqC.
% sort(v)
% 
% [Command Window outputs]
% | ans =
% |     1.6209    1.9528    3.2386    3.3837
% | ans =
% |     1.6209    1.9528    3.2386    3.3837
% 
% While we're here, let's use the "matrix square root" to make some random
% draws and then examine the covariance.
% 
% d = sqC * randn(4, 100000);
% C_empirical = cov(d.')
% C
% 
% [Command Window outputs]
% | C_empirical =
% |     2.1932   -0.1192   -0.1265   -0.5776
% |    -0.1192    2.6547   -0.7807   -0.0465
% |    -0.1265   -0.7807    2.3431    0.0947
% |    -0.5776   -0.0465    0.0947    3.0091
% | C =
% |     2.1788   -0.1236   -0.1120   -0.5674
% |    -0.1236    2.6326   -0.7939   -0.0444
% |    -0.1120   -0.7939    2.3581    0.0911
% |    -0.5674   -0.0444    0.0911    3.0265
% 
% See also: covdraw

% Copyright 2017 An Uncommon Lab
    
    % Defaults
    if nargin < 1 || isempty(n),       n       = 10;    end;
    if nargin < 2 || isempty(eig_min), eig_min = 0;     end;
    if nargin < 3 || isempty(eig_max), eig_max = 1;     end;
    if nargin < 4 || isempty(type),    type    = 'cov'; end;

    % Create an orthonormal matrix. We'll use QR to make life easy.
    [C, ~] = qr(randn(n));
    
    % Draw some eigenvalues.
    v = eig_min + (eig_max - eig_min) .* rand(1, n);
    
    % Construct the covariance matrix as a translation of the eigenvalues.
    switch lower(type)
        case 'cov'
            C = C * diag(v) * C.';
        case 'sqrt'
            C = C * diag(sqrt(v));
        otherwise
            error('Unknown covariance type: %s.', type);
    end
    
    % If there's no output, show some debug stuff.
    if nargout == 0
        
        % Draw some samples from a multivariate normal distribution with
        % this covariance matrix.
        switch lower(type)
            case 'cov'
                x = mnddraw(C, 1e6);
            case 'sqrt'
                x = C * randn(n, 1e6);
                C = C * C.';
        end
        
        % Determine the covariance from the samples.
        C_empirical = cov(x.');
        
        % Show the results.
        figure();
        subplot(2, 1, 1);
        imagesc(C);
        title('Random Covariance Matrix');
        colorbar();
        subplot(2, 1, 2);
        imagesc(C_empirical);
        colorbar();
        title('Empirical Covariance Matrix');
        
        fprintf('Drawn eigenvalues: \n');
        sort(v)
        fprintf('Calculated eigenvalues: \n');
        sort(eig(C).')
        fprintf('Empirical eigenvalues: \n');
        sort(eig(C_empirical).')
        
    end % example
    
end % randcov

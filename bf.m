function [x_hat_k, X_k, w_k, w_kkm1] = bf( ...
    t_km1, t_k, X_km1, w_km1, u_km1, z_k,  ...
    f, h, draw_proc_noise, prob_obs_error, ...
    h_tune, resample)
                              
% bf
% 
% Run one update of a basic boostrap filter, a type of particle filter, 
% which can estimate well enough on many types of problems.
%
%    [x_k, X_k, w_k] = bf(t_km1, t_k, X_km1, w_km1, u_km1, z_k, ...
%                         f, h, draw_proc_noise, prob_obs_error, ...
%                         h_tune);
% 
% Inputs:
%
% t_km1            Time corresponding to the input states (sample k-1)
% t_k              Time at the current measurement (sample k)
% X_km1            Matrix of particles at sample k-1, |nx-by-np|
% w_km1            Matrix of weights at sample k-1, |1-by-np|
% u_km1            Input vector at sample k-1
% z_k              Measurement at sample k
% f                Propagation function with an interface defined as:
%
%                     x_k = fcn([t_km1 t_k], x_km1, u_km1, q_km1)
%
%                  where |x_km1| is the prior state of a single particle
%                  and |q_km1| is the process noise acting on the
%                  particle, yielding an updated state, |x_k|, for the
%                  particle, and |t_k| is the time of current sample
% h                Measurement function with an interface defined as:
%
%                     z_k = fcn(t_k, x_k, u_km1)
%
%                  where |x_k| is the state of a particle at sample k and
%                  |z_k| is the corresponding measurement
% draw_proc_noise  A function to draw process noise given a time and
%                  corresponding particle state
%
%                     q_km1 = fcn(t_km1, x_km1)
% 
% prob_obs_error   A function to return the probability density of the
%                  measurement error -- the difference between the
%                  actual measurement and a perfect measurement for the
%                  true state, with the interface:
%
%                     p_k = fcn(t_k, dz_k, x_k)
%
%                  Only relative scores matter, so constant coefficients
%                  can be dropped. For instance, for a measurement model
%                  with normally distrubted errors with covariance
%                  matrix, |R|, one can use:
%
%                     p_k = @(t, dz, x) exp(-0.5 * dz.' * inv(R) * dz)
%
%                  (noting that the preceding constant
%                  |(2*pi)^(-k/2) / sqrt(det(R))| has been dropped since
%                  it won't affect the results).
% h_tune           An optional regulation tuning parameter, defaulting
%                  to:
%
%                     h_tune = (4/(N * (nx + 2)))^(1/(nx + 4));
%
%                  which is optimal for a Gaussian error distribution.
%                  This parameter scales the spread of particles when
%                  they are redrawn about the mean during regularization.
% resample         Set to true to resample and regularize (true by default)
% 
% Outputs:
%
% x_hat_k       Mean state estimate for sample k
% X_k           Set of particles for sample k
% w_k           Particle weights at sample k
% w_kkm1        Updated weights, prior to resampling (for analysis only)
% 
% See <https://www.anuncommonlab.com/doc/starkf/bf.html> for more.

% Copyright 2016 An Uncommon Lab

    % Dimensions
    [nx, N] = size(X_km1);
    nz = size(z_k, 1);
    
    % If no tuning parameter is provided, use the optimal regularization 
    % parameter for a Gaussian kernal.
    if nargin < 11 || isempty(h_tune)
        h_tune = (4/(N * (nx + 2)))^(1/(nx + 4));
    end

    % Resample and regularize by default.
    if nargin < 12, resample = true; end;

    % Propagate.
    X_kkm1 = zeros(nx, N);
    Z_kkm1 = zeros(nz, N);
    w_kkm1 = zeros(1,  N);
    for k = 1:N
        
        % Predict the state.
        X_kkm1(:, k) = f(t_km1, t_k, X_km1(:,k), u_km1, ...
                         draw_proc_noise(t_km1, X_km1(:,k)));
                     
        % Predict the measurement.
        Z_kkm1(:, k) = h(t_k, X_kkm1(:,k), u_km1);
        
        % Update the weight based on the difference between the true and
        % predicted measurement.
        w_kkm1(k) = w_km1(k) * prob_obs_error(t_k, ...
                                              z_k - Z_kkm1(:, k), ...
                                              X_kkm1(:, k));

    end
    
    % Normalize the weights.
    w_kkm1 = w_kkm1 ./ sum(w_kkm1);
    
    % Calculate the updated state. This now includes information from the
    % observation since we updated the weights with the probability of each
    % observation error.
    x_hat_k = sum(bsxfun(@times, w_kkm1, X_kkm1), 2);

    % If we're resampling and regularizing...
    if resample

        % Randomly pick particles from X_kkm1 according to their weights.
        % The weights act as a probability mass function.
        indices = pmfdraw(w_kkm1);
        X_k     = X_kkm1(:, indices);
        
        % By virtue of randomly selecting particles, all of the current
        % particles now have the same weight of 1/N.
        w_k = repmat(1/N, 1, N);

        % Regularize by creating the sample covariance matrix and drawing
        % random perturbations from it.
        X_twiddle_k = bsxfun(@minus, X_k, sum(X_k, 2)/N);
        Sigma       = 1/(N-1) * (X_twiddle_k * X_twiddle_k.');
        
%         % Alternately, we could have used the weights instead of the
%         % newly selected particles. The expected result is the same; the
%         % actual result may be slightly different.
%         X_twiddle_k = bsxfun(@minus, X_kkm1, x_hat_k);
%         X_twiddle_k = bsxfun(@times, X_twiddle_k, sqrt(w_kkm1));
%         Sigma       = X_twiddle_k * X_twiddle_k.';
        
        % Multiply those perturbations by the tuning parameter and add them
        % to the particles.
        [u, s] = svd(Sigma);
        X_k    = X_k + h_tune * u * sqrt(s) * randn(nx, N);
        
    % Otherwise, use the updated values.
    else
        X_k = X_kkm1;
        w_k = w_kkm1;
    end

end % bf

% Draw n indices randomly from a probability mass function defined by the
% vector p.
function indices = pmfdraw(p, n)

    % Default number of draws is the number of elements in p.
    if nargin < 2, n = length(p); end;

    % Create the discrete cumulative distribution function.
    cdf = cumsum(p);
    
    % If it's not actually a normalized input, normalize it.
    if cdf(end) ~= 1
        cdf = cdf ./ cdf(end);
    end
    
    % Make a bunch of normal, uniform random draws and sort them.
    draws = sort(rand(1, n));
    
    % Find the CDF "bin" containing each draw. This will be the appropriate
    % index of the CDF.
    np      = length(p);
    indices = zeros(1, n);  % Output randomly drawn indices
    ci      = 1;            % Index of PMF bin
    for di = 1:n
        while cdf(ci) < draws(di) && ci < np
            ci = ci + 1;
        end
        indices(di) = ci;
    end
    
end % pmfdraw

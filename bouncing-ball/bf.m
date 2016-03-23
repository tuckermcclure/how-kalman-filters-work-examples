function [x_hat_k, X_k, w_k, w_kkm1] = bf(t_km1, t_k, X_km1, w_km1, u_km1, z_k, ...
                                  f, h, draw_proc_noise, prob_obs_error,...
                                  h_tune, resample)
                              
%#autodoc

% Copyright 2015 An Uncommon Lab

%#codegen

    % Dimensions
    nx = size(X_km1, 1);
    nz = size(z_k,   1);
    N  = size(X_km1, 2);
    
    % Optimal regularization parameter for Gaussian kernal
    if nargin < 11 || isempty(h_tune)
        h_tune = (4/(N * (nx + 2)))^(1/(nx + 4));
    end

    % Predict.
    X_kkm1 = zeros(nx, N);
    Z_kkm1 = zeros(nz, N);
    w_k    = zeros(1,  N);
%     dz = zeros(nz, N);
    for k = 1:N
        X_kkm1(:, k) = f(t_km1, t_k, X_km1(:,k), u_km1, ...
                         draw_proc_noise(t_km1, X_km1(:,k)));
        Z_kkm1(:, k) = h(t_k, X_kkm1(:,k), u_km1);
        w_k(k)       = w_km1(k) * prob_obs_error(t_k, ...
                                                 z_k - Z_kkm1(:, k), ...
                                                 X_kkm1(:, k));
%         dz(:, k) = z_k - Z_kkm1(:, k);
    end
    w_k = w_k ./ sum(w_k);
    w_kkm1 = w_k;
    
    % Calculate the updated state. This now includes information from the
    % observation since we updated the weights with the probability of each
    % observation error.
    x_hat_k = sum(bsxfun(@times, w_k, X_kkm1), 2);

    if resample

        % Update weights and resample. By virtue of resampling, all points 
        % are now equally likely, so all w_k == 1/N.
%         w_k = sqrt(w_k);
%         w_k = w_k / sum(w_k);
        indices = pmfdraw(w_k); % Discretely sample X_kkm1 according to w_k.
        X_k     = X_kkm1(:, indices);
        w_k(:)  = 1/N;

        % Regularize.
        X_twiddle_k = bsxfun(@minus, X_k, sum(X_k, 2)/N);
        Sigma       = 1/(N-1) * (X_twiddle_k * X_twiddle_k.');
        X_k         = X_k + h_tune * sqrtpsdm(Sigma) * randn(nx, N);
        
    else
        X_k = X_kkm1;
    end

end % bf

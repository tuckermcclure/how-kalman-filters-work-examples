function [x, P, x_kkm1, P_kkm1, P_xz, P_zz, K] = ukfan(t_km1, t_k, x, P, u_km1, z_k, ...
                        f, h, Q_km1, R_k, ...
                        alpha, beta, kappa, varargin)

% ukfan
%
% Run one update of a basic unscented (sigma point) Kalman filter where the
% process and measurement noise are additive. This runs faster than a
% traditional UKF and, because the propagation and observation functions do
% not require inputs for the process and measurement noise, is a drop-in
% replacement for an extended Kalman filter.
%
%   [x_k, P_k] = ukfan( ...
%       t_km1, t_k, x_km1, P_km1, u_km1, z_k, ...
%       f, h, Q_km1, R_k, ...
%       alpha, beta, kappa, ...
%       (etc.))
%
% 
% Inputs:
%
% t_km1    Time at sample k-1
% dt       Time step from sample k-1 to k
% x_km1    State estimate at sample k-1
% P_km1    Estimate covariance at sample k-1
% u_km1    Input vector at sample k-1
% z_k      Measurement at sample k
% f        Propagation function with interface:
%
%            x_k = f(t_km1, t_k, x_km1, u_km1, q_km1)
%
%          where |q_km1| is the process noise at sample k-1
%
% h        Measurement function with interface:
%
%            z_k = h(t_k, x_k, u_km1, r_k)
%
%          where |r_k| is the measurement noise at sample k
%
% Q_km1    Pocess noise covariance at sample k-1
% R_k      Measurement noise covariance at sample k
% alpha    Optional turning parameter, often 0.001
% beta     Optional tuning parameter, with 2 being optimal for Gaussian
%          estimation error
% kappa    Optional tuning parameter, often 3 - |nx|, where |nx| is the
%          dimension of the state
% (etc.)   Additional arguments to be passed to f and h after their
%          normal arguments
% 
% Outputs:
%
% x        Upated estimate at sample k
% P        Updated estimate covariance at sample k
% 
% Example:
% 
% We can quickly create a simulation for discrete, dynamic system, generate
% noisy measurements of the system over time, and pass these to an
% unscented Kalman filter.
% 
% First, define the discrete system.
% 
% rng(1);
% dt    = 0.1;                                  % Time step
% F_km1 = expm([0 1; -1 0]*dt);                 % State transition matrix
% G_km1 = [0.5*dt^2; dt];                       % Process-noise-to-state map
% Q_km1 = G_km1 * 0.5^2 * G_km1.';              % Process noise variance
% R_k   = 0.1;                                  % Meas. noise covariance
% 
% Make propagation and observation functions. These are just linear for
% this example, but |ukf| is meant for nonlinear functions.
% 
% f = @(t_km1, t_k, x_km1, u_km1) F_km1 * x_km1;
% h = @(t_k, x_k, u_k) x_k(1);
% 
% Now, we'll define the simulation's time step and initial conditions. Note
% that we define the initial estimate and set the truth as a small error
% from the estimate (using the covariance).
% 
% n       = 100;                     % Number of samples to simulate
% x_hat_0 = [1; 0];                  % Initial estimate
% P       = diag([0.5 1].^2);        % Initial estimate covariance
% x_0     = x_hat_0 + mnddraw(P, 1); % Initial true state
% 
% Now we'll just write a loop for the discrete simulation.
% 
% % Storage for time histories
% x     = [x_0, zeros(2, n-1)];                         % True state
% x_hat = [x_hat_0, zeros(2, n-1)];                     % Estimate
% z     = [x_0(1) + mnddraw(R_k, 1), zeros(1, n-1)];    % Measurement
% 
% % Simulate each sample over time.
% for k = 2:n
% 
%     % Propagate the true state.
%     x(:, k) = F_km1 * x(:, k-1) + mnddraw(Q_km1, 1);
%     
%     % Create the real measurement at sample k.
%     z(:, k) = x(1, k) + mnddraw(R_k, 1);
% 
%     % Run the Kalman correction.
%     [x_hat(:,k), P] = ukfan((k-1)*dt, dt, x_hat(:,k-1), P, [], z(:,k), ...
%                             f, h, Q_km1, R_k);
% 
% end
% 
% Plot the results.
% 
% figure(1);
% clf();
% t = 0:dt:(n-1)*dt;
% plot(t, x, ...
%      t, z, '.', ...
%      t, x_hat, '--');
% legend('True x1', 'True x2', 'Meas.', 'Est. x1', 'Est. x2');
% xlabel('Time');
% 
% Note how similar this example is to the example of |ukf|.
% 
% 
% 
% Reference
%
% Wan, Eric A. and Rudoph van der Merwe. "The Unscented Kalman Filter."
% _Kalman Filtering and Neural Networks._ Ed. Simon Haykin. New York: John
% Wiley & Sons, Inc., 2001. Print. Pages 221-276.
% 
% See also: srukf, ukf

                    
% Copyright 2015 An Uncommon Lab

%#codegen

    if nargin < 11
        alpha = 0.001;
        beta  = 2;
        kappa = 3 - numel(x);
    end

    % Dimensions
    nx = length(x);      % Number of states
    nz = size(R_k, 1);   % Number of observations
    L  = nx;             % Dimension of augmented state
    ns = 2 * L + 1;      % Number of sigma points
    
    % Weights
    lambda = alpha^2 * (L + kappa) - L;
    gamma  = sqrt(L + lambda);
    W_m_0  = lambda / (L + lambda);
    W_c_0  = W_m_0 + 1 - alpha^2 + beta;
    W_i    = 1/(2*(L + lambda));
    
    % Create the sigma points.
    gamma_sqrt_P = gamma * sqrtpsdm(P);
    X_km1 = [x, ...
             repmat(x, 1, L) + gamma_sqrt_P, ...
             repmat(x, 1, L) - gamma_sqrt_P];

    % Predict sigma points and measurements.
    X_star_kkm1 = zeros(nx, ns);
    for k = 1:ns
        X_star_kkm1(:, k) = f(t_km1, t_k, X_km1(:, k), u_km1, varargin{:});
    end
    
    % Expected prediction and measurement
    x_kkm1 = W_m_0 * X_star_kkm1(:,1) + W_i * sum(X_star_kkm1(:,2:end), 2);
    
    % We don't need X_star_kkm1 any longer, so remove the expectations in
    % order to calculate the covariance.
    X_star_kkm1 = bsxfun(@minus, X_star_kkm1, x_kkm1);
    
    % Calculate covariance of the prediction.
    P_kkm1 =   (W_c_0 * X_star_kkm1(:, 1)) * X_star_kkm1(:, 1).' ...
             + W_i * (X_star_kkm1(:, 2:end) * X_star_kkm1(:, 2:end).') ...
             + Q_km1;
          
    % Make the augmented sigma points. Note that L grows, and so gains
    % change.
    L            = 2 * L;
    ns           = 2 * L + 1;
    lambda       = alpha^2 * (L + kappa) - L;
    gamma_old    = gamma; % We'll need this to scale the old sigma points.
    gamma        = sqrt(L + lambda);
    gamma_sqrt_Q = gamma * sqrtpsdm(Q_km1);

    % Calculate the new weights for the augmented sigma points.
    W_m_0  = lambda / (L + lambda);
    W_c_0  = W_m_0 + 1 - alpha^2 + beta;
    W_i    = 1/(2*(L + lambda));

	% Scale the old sigma points' displacements from x_kkm1 for the new 
    % gains.
    X_star_kkm1 = X_star_kkm1 * gamma/gamma_old;
    
	% Create the new sigma points which now include the process noise.
    X_kkm1 = [repmat(x_kkm1, 1, 2*nx+1) + X_star_kkm1, ...
              repmat(x_kkm1 + X_star_kkm1(:, 1), 1, nx) + gamma_sqrt_Q, ...
              repmat(x_kkm1 + X_star_kkm1(:, 1), 1, nx) - gamma_sqrt_Q];

    % Predict sigma points and measurements.
    Z_kkm1 = zeros(nz, ns);
    for k = 1:ns
        Z_kkm1(:, k) = h(t_k, X_kkm1(:, k), u_km1, varargin{:});
    end
    
    % Calculate the expectation and remove it from the set.
    z_kkm1 = W_m_0 * Z_kkm1(:, 1) + W_i * sum(Z_kkm1(:, 2:end), 2);
    Z_kkm1 = bsxfun(@minus, Z_kkm1, z_kkm1);
    X_kkm1 = bsxfun(@minus, X_kkm1, x_kkm1);
         
    % Covariance of predicted observation
    P_zz =   (W_c_0 * Z_kkm1(:, 1)) * Z_kkm1(:, 1).' ...
           + W_i * (Z_kkm1(:, 2:end) * Z_kkm1(:, 2:end).') ...
           + R_k;
    
    % Covariance of predicted observation and predicted state
    P_xz =   (W_c_0 * X_kkm1(:, 1)) * Z_kkm1(:, 1).' ...
           + W_i * (X_kkm1(:, 2:end) * Z_kkm1(:, 2:end).');
    
    % Kalman gain
    K = P_xz / P_zz;
    
    % Correct the state.
    if isempty(z_k)
        
        % z_kkm1 is the innovation vector, and so K has opposite sign of 
        % the intended gain due to the P_xz calculation. This won't matter
        % for the covariance below (the sign is squared).
        x = x_kkm1 - K * z_kkm1;
        
    else
        x = x_kkm1 + K * (z_k - z_kkm1);
    end
    
    % Correct the covariance.
    P = P_kkm1 - K * P_zz * K.';
    
end % ukfan

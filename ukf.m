function [x, P, x_kkm1, P_kkm1, P_xz, P_zz, K] = ukf( ...
    t_km1, t_k, x, P, u_km1, z_k, ...
    f, h, Q_km1, R_k, ...
    alpha, beta, kappa, ...
    varargin)

% ukf
%
% Run one update of a basic unscented (sigma point) Kalman filter. The
% advantage of a UKF over a traditional extended Kalman filter is that the
% UKF can be more accurate when the propagation/measurement functions are
% not zero-mean wrt the error (that is, they are significantly nonlinear in
% terms of their errors). Further, since a UKF doesn't require Jacobians,
% it's often easier to use a UKF than an EKF. The disadvantage is increased
% runtime. Though a UKF operates "on the order" of an EKF, in
% implementation it often takes several times longer.
%
%   [x_k, P_k] = ukf(t_km1, t_k, x_km1, P_km1, u_km1, z_k, ...
%                    f, h, Q_km1, R_k, ...
%                    alpha, beta, kappa, ...
%                    (etc.))
% 
% Inputs:
%
% t_km1    Time at sample k-1
% t_k      Time at sample k
% x_km1    State estimate at sample k-1
% P_km1    Estimate covariance at sample k-1
% u_km1    Input vector at sample k-1
% z_k      Measurement at sample k
% f        Propagation function with interface:
%
%            x_k = f([t_km1, t_k], x_km1, u_km1, q_km1)
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
% alpha    Optional turning parameter (e.g., 0.001)
% beta     Optional tuning parameter, with 2 being optimal for Gaussian
%          estimation error
% kappa    Optional tuning parameter, often 3 - |nx|, where |nx| is the
%          dimension of the state
% (etc.)   Additional arguments to be passed to |f| and |h| after their
%          normal arguments
% 
% Outputs:
%
% x        Upated estimate at sample k
% P        Updated estimate covariance at sample k
% 
% Reference:
%
% Wan, Eric A. and Rudoph van der Merwe. "The Unscented Kalman Filter."
% _Kalman Filtering and Neural Networks._ Ed. Simon Haykin. New York: John
% Wiley & Sons, Inc., 2001. Print. Pages 221-276.
% 
% See <https://www.anuncommonlab.com/doc/starkf/ukf.html> for more.

% Copyright 2016 An Uncommon Lab
          
	% Defaults
    if nargin < 11 || isempty(alpha), alpha = 0.001;        end;
    if nargin < 12 || isempty(beta),  beta  = 2;            end;
    if nargin < 13 || isempty(kappa), kappa = 3 - numel(x); end;
                        
    % Dimensions
    nx = length(x);         % Number of states
    nq = size(Q_km1, 1);    % Number of noise states
    nz = size(R_k, 1);      % Number of observations
    L  = nx + nq + nz;      % Dimension of augmented state
    ns = 2 * L + 1;         % Number of sigma points
    
    % Weights
    lambda = alpha ^2 * (L + kappa) - L;
    gamma  = sqrt(L + lambda);
    W_m_0  = lambda / (L + lambda);
    W_c_0  = W_m_0 + 1 - alpha^2 + beta;
    W_i    = 1/(2*(L + lambda));

    % Create the augmented system.
    P_a     = blkdiag(P, Q_km1, R_k);
    x_a_km1 = [x; zeros(nq, 1); zeros(nz, 1)];
    
    % Create the sigma points. (We use svd instead of chol in case there's
    % a 0 eigenvalue.)
    [u, s] = svd(P_a);
    gamma_sqrt_P_a = gamma * u * sqrt(s);
    X_a_km1 = [x_a_km1, ...
               repmat(x_a_km1, 1, L) + gamma_sqrt_P_a, ...
               repmat(x_a_km1, 1, L) - gamma_sqrt_P_a];

    % Predict sigma points and measurements.
    Z_kkm1   = zeros(nz, ns);
    X_x_kkm1 = zeros(nx, ns);
    for sp = 1:ns
        X_x_kkm1(:, sp) = f(t_km1, t_k, ...              % Times
                            X_a_km1(1:nx, sp), ...       % State k-1
                            u_km1, ...                   % Input
                            X_a_km1(nx+1:nx+nq, sp), ... % Proc. noise
                            varargin{:});                % Et al.
        Z_kkm1(:, sp)   = h(t_k, ...                     % Time k
                            X_x_kkm1(:, sp), ...              % State k|k-1
                            u_km1, ...                        % Input
                            X_a_km1(nx+nq+1:nx+nq+nz, sp), ...% Meas. noise
                            varargin{:});                     % Et al.
    end
    
    % Expected prediction and measurement
    x_kkm1 = W_m_0 * X_x_kkm1(:, 1) + W_i * sum(X_x_kkm1(:, 2:end), 2);
    z_kkm1 = W_m_0 * Z_kkm1(:, 1)   + W_i * sum(Z_kkm1(:, 2:end), 2);

    % Remove expectations from X_x_kkm1 and Z_kkm1.
    X_x_kkm1 = bsxfun(@minus, X_x_kkm1, x_kkm1);
    Z_kkm1   = bsxfun(@minus, Z_kkm1, z_kkm1);
    
    % Calculate covariance of the prediction.
    P_kkm1 =   (W_c_0 * X_x_kkm1(:, 1)) * X_x_kkm1(:, 1).' ...
             + W_i * (X_x_kkm1(:, 2:end) * X_x_kkm1(:, 2:end).');
    
    % Covariance of predicted observation
    P_zz =   (W_c_0 * Z_kkm1(:, 1)) * Z_kkm1(:, 1).' ...
           + W_i * (Z_kkm1(:, 2:end) * Z_kkm1(:, 2:end).');
    
    % Covariance of predicted observation and predicted state
    P_xz =   (W_c_0 * X_x_kkm1(:, 1)) * Z_kkm1(:, 1).' ...
           + W_i * (X_x_kkm1(:, 2:end) * Z_kkm1(:, 2:end).');
       
    % Kalman gain
    K = P_xz / P_zz;
    
    % Correct the state.
    x = x_kkm1 + K * (z_k - z_kkm1);
    
    % Correct the covariance.
    P = P_kkm1 - K * P_zz * K.';
    
end % ukf

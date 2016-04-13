function [x_hat, P] = ekf(t_km1, t_k, x_hat, P, u, z, ...
                          f, h, F_fcn, H_fcn, Q, R)

% ekf
% 
% Runs a single iteration of an extended Kalman filter.
%
% Inputs:
%
% t_km1    Time at sample k-1
% t_k      Time at sample k
% x_hat    Estimate at sample k-1 (nx-by-1)
% P        Estimate covariance matrix at sample k-1 (nx-by-nx)
% u        Input vector (nu-by-1)
% z        Measurement (nz-by-1)
% f        Handle of the propagation function, mapping a state at time
%          t_km1 to t_k, with the following interface:
%
%             x_k = f(t_km1, t_k, x_km1, u);
%
% h        Observation function, mapping a state to a predicted
%          measurement, with the following interface:
% 
%             z_hat_k = h(t_k, x_k, u);
%
% F_fcn    Function producing the Jacobian of the propagation function from
%          time t_km1 to t_k:
%
%             F = F_fcn(t_km1, t_k, x_km1, u);
%
% H_fcn    Function producing the Jacobian of the observation function:
%
%             H = H_fcn(t_k, x_k, u);
%
% Q        Process noise covariance matrix (nx-by-nx)
% R        Measurement noise covariance matrix (nz-by-nz)
% 
% Outputs:
%
% x_hat    Estimate at sample k
% P        Estimate covariance matrix at sample k
% 

% Copyright 2016 An Uncommon Lab

    % Calculate the Jacobian using the state at k-1.
    F = F_fcn(t_km1, t_k, x_hat, u);
    
    % Propagate the estimate.
    x_hat = f(t_km1, t_k, x_hat, u);
    
    % Propagate the covariance.
    P = F * P * F.' + Q;
    
    % Predict the measurement using the predicted state estimate.
    z_hat = h(t_k, x_hat, u);
    
    % Innovation vector
    y = z - z_hat;
    
    % Calculate the observation Jacobian using the predicted state
    % estimate.
    H = H_fcn(t_k, x_hat, u);
    
    % Calculate the state-innovation covariance and innovation covariance.
    P_xy = P * H.';
    P_yy = H * P * H.' + R;
    
    % Calculate the Kalman gain.
    K = P_xy / P_yy;
    
    % Correct the estimate.
    x_hat = x_hat + K * y;
    
    % Correct the covariance using Joseph form for stability. This is the
    % same as P = P - K * H * P, but presents less of a problem in the
    % presense of floating point roundoff.
    A = eye(length(x_hat)) - K * H;
    P = A * P * A.' + K * R * K.';

end % ekf

function [x_hat, P] = lkf(x_hat, P, u, z, F, F_u, H, H_u, Q, R)

% lkf
% 
% Runs a single iteration of a linear Kalman filter.
%
% Inputs:
%
% x_hat    Estimate at sample k-1 (nx-by-1)
% P        Estimate covariance matrix at sample k-1 (nx-by-nx)
% u        Input vector (nu-by-1)
% z        Measurement (nz-by-1)
% F        State transition matrix (nx-by-nx)
% F_u      Map from input vector to state (nx-by-nu)
% H        Observation matrix (nz-by-nx)
% H_u      Map from input vector to measurement (nz-by-nu)
% Q        Process noise covariance matrix (nx-by-nx)
% R        Measurement noise covariance matrix (nz-by-nz)
% 
% Outputs:
%
% x_hat    Estimate at sample k
% P        Estimate covariance matrix at sample k
% 

% Copyright 2016 An Uncommon Lab

    % Propagate the estimate.
    x_hat = F * x_hat + F_u * u;
    
    % Propagate the covariance.
    P = F * P * F.' + Q;
    
    % Predict the measurement using the predicted state estimate.
    z_hat = H * x_hat + H_u * u;
    
    % Innovation vector
    y = z - z_hat;
    
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

end % lkf

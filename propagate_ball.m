function [x, xt, tt] = propagate_ball(t_km1, t_k, x, varargin)

% propagate_ball
%
% Propagate the ball forward in time by dt, bouncing if necessary.
%
% Inputs:
%
% t_km1   Start time of simulation [s]
% t_k     Stop time of simulation [s]
% x       State of ball ([x; y; x_dot; y_dot]) at t_km1
%
% Outputs:
%
% x       Propagated state at t_k
% xt      Trajectory from t_km1 to t_k (for analysis only)
% tt      Trajectory sample times (for analysis only)

% Copyright 2016 An Uncommon Lab

    % Constants
    g = 9.81;
    
    % Start the trajectory.
    xt = x;
    tt = t_km1;

    % For each minor time step...
    t = t_km1;
    while t < t_k - eps

        % Calculate the time step.
        dt = min(0.1, t_k - t);
        t = t + dt;
    
        % Calculate the time until the next bounce. Since the parabola
        % always faces downwards, the bounce will be the larger of the two
        % solutions to the quadratic equation:
        % 
        %   -0.5 * g * dT^2 + x(4) * dT + x(2) = 0;
        % 
        % However, if a particle ends up below ground and going down, there
        % will be no positive solution.
        if x(4)^2 + 2 * g * x(2) < 0
            dT = 0;
        else
            dT = (x(4) + sqrt(x(4)^2 + 2 * g * x(2))) / g;
        end
        
        % We never go back in time to bounce.
        dT = max(0, dT);
        
        % If we bounce on this time step...
        if dT < dt

            % x(2) = x(2) + x(4) * dT - 0.5 * 9.81 * dT^2; % Goes to 0.
            x(4) = -0.95 * (x(4) - g * dT);           % Bounce after dT.
            xi   = [x(1) + x(3) * dT; 0; x(3); x(4)]; % Bounce state
            ti   = t - dt + dT;                       % Time of bounce
            dT   = dt - dT; % Time from bounce to the end of the time step
            x(2) = 0 + x(4) * dT - 0.5 * 9.81 * dT^2; % Bounce up from 0.
            x(4) = x(4) - 9.81 * dT; % Update velocity for remainder of step.

        % Otherwise, no bounce.
        else

            x(2) = x(2) + x(4) * dt - 0.5 * g * dt^2;
            x(4) = x(4) - 9.81 * dt;
            xi   = [];
            ti   = [];

        end
        
        % We don't check for a second bounce on a time step. We could, but
        % it wouldn't matter for this example.
        
        % Updating the x direction is trivial and the same regardless of
        % bounce.
        x(1) = x(1) + x(3) * dt;
        % x(3) = x(3); % The x velocity never changes.
        
        % Create the trajectory if necessary.
        if nargout >= 2
            xt = [xt, xi, x]; %#ok<AGROW>
            tt = [tt, ti, t]; %#ok<AGROW>
        end
        
    end % each minor time step
    
end % propagate_ball

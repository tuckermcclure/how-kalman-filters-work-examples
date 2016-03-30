function [x, xt, tt] = propagate(t_km1, t_k, x, varargin)

% propagate
%
% Propagate the ball forward in time by dt, bouncing if necessary.
%
% Inputs:
%
% x     State of ball ([x; y; x_dot; y_dot])
% dt    Time to propagate
%
% Outputs:
%
% x     Propagated state

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
    
        % Calculate the time until bounce.
        % -0.5 * 9.81 * dT^2 + x(4) * dT + x(2) = 0;
        if x(4)^2 + 2 * g * x(2) < 0
            dT = 0;
        else
            dT = (x(4) + sqrt(x(4)^2 + 2 * g * x(2))) / g;
        end
        
        % If we bounce on this time step...
        if dT < dt

            % x(2) = x(2) + x(4) * dT - 0.5 * 9.81 * dT^2; % Goes to 0.
            x(4) = -0.95 * (x(4) - g * dT); % Bounce after dT.
            xi   = [x(1) + x(3) * dT; 0; x(3); x(4)];
            ti   = t - dt + dT;
            dT   = dt - dT; % Time from bounce to the end fo the time step
            x(2) = 0 + x(4) * dT - 0.5 * 9.81 * dT^2; % Bounce up from floor.
            x(4) = x(4) - 9.81 * dT; % Update velocity for remainder of step.

        % Otherwise, no bounce.
        else

            x(2) = x(2) + x(4) * dt - 0.5 * g * dt^2;
            x(4) = x(4) - 9.81 * dt;
            xi   = [];
            ti   = [];

        end

        % Updating the x direction is trivial and the same regardless of
        % bounce.
        x(1) = x(1) + x(3) * dt;
        % x(3) = x(3);
        
        % Create the trajectory if necessary.
        if nargout >= 2
            xt = [xt, xi, x]; %#ok<AGROW>
            tt = [tt, ti, t]; %#ok<AGROW>
        end
        
    end % each minor time step
    
    if any(imag(x))
        disp('hi');
    end
    
end % propagate

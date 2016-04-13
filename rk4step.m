function x_k = rk4step(f, t_km1, x_km1, dt, varargin)

% rk4step
%
% Runs a single step of the common fourth order Runge-Kutta numerical
% integration method.
%
% Inputs:
%
% f      Function handle taking time, state, and any additional arguments
%        and returning state time rate of change:
%
%           x_dot = f(t, x, ...)
% 
% t_km1  Time at sample k-1
% x_km1  State at sample k-1
% dt     Time step from sample k-1 to k
%
% Outputs:
%
% x_k    State at sample k
%

% Copyright 2016 An Uncommon Lab

%#codegen

    % Use RK4 to propagate.
    d1    = f(t_km1,          x_km1,             varargin{:});
    d2    = f(t_km1 + 0.5*dt, x_km1 + 0.5*dt*d1, varargin{:});
    d3    = f(t_km1 + 0.5*dt, x_km1 + 0.5*dt*d2, varargin{:});
    d4    = f(t_km1 +     dt, x_km1 +     dt*d3, varargin{:});
    x_k = x_km1 + dt/6*(d1 + 2*d2 + 2*d3 + d4);

end % rk4step

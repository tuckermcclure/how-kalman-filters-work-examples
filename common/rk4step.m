function x_kp1 = rk4step(f, t_k, x_k, dt, varargin)

% rk4step
%
% Runs a single step of the fourth order Runge-Kutta numerical integration
% method.
%
% Inputs:
%
% f     Function handle taking time, state, and any additional arguments
%       and returning state time rateof change:
%
%          x_dot = f(t, x, ...)
% 
% t_k   Time at sample k
% x_k   State at sample k
% dt    Time step from sample k to k+1
%
% Outputs:
%
% x_kp1 State at sample k+1
%
% Copyright 2016 An Uncommon Lab

    % Use RK4 to propagate.
    d1    = f(t_k,          x_k,             varargin{:});
    d2    = f(t_k + 0.5*dt, x_k + 0.5*dt*d1, varargin{:});
    d3    = f(t_k + 0.5*dt, x_k + 0.5*dt*d2, varargin{:});
    d4    = f(t_k +     dt, x_k +     dt*d3, varargin{:});
    x_kp1 = x_k + dt/6*(d1 + 2*d2 + 2*d3 + d4);

end % rk4step

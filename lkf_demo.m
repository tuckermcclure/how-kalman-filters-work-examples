%% Linear Kalman Filter Demonstration
%
% This script demonstrates using a linear Kalman filter to determine
% the position and velocity of a falling package for the article "How
% Kalman Filters Work".
%
% <http://www.anuncommonlab.com/articles/how-kalman-filters-work/>
%
% Copyright 2016 Tucker McClure @ An Uncommon Lab

%% System definition and simulation
%
% We'll define the system, set up the filter, and simulate the results over
% 10s.

% Set the random number generator seed so the results are the same every
% time we run the script. (Comment out this line to see different results
% every time.)
rng(8);

% Define the true system.
dt = 0.1;        % Time step [s]
cd = 4.4;        % Coefficient of drag, so that F_d = -cd * norm(v) * v;
m  = 1;          % Mass [kg]
g  = [0; -9.81]; % Gravity [m/s^2]

% Define the continuous-time dynamics, x_dot = f(x).
fc = @(t, x, q) [x(3:4); g - cd/m * norm(x(3:4)) * x(3:4) + q];

% Define a function to update the state by 1 time step with the given
% process noise. When the filter uses this, the process noise will
% naturally be [0; 0].
f = @(x, q) rk4step(fc, 0, x, dt, q);

% Define the process and measurement noise.
Q  = 0.5^2 * eye(2);  % Acceleration noise, 0.5m/s^2 std. dev.
R  = 0.5^2 * eye(2);  % 0.5m standard deviation

% Initial conditions
x0  = [-1; 8; 1; 0];            % True state [m; m; m/s; m/s]
z0  = x0(1:2) + covdraw(R);     % Initial measurement [m; m]

% The initial estimate of position will come directly from the first
% measurement, so the corresponding part of the initial covariance will
% have the same covariance as the measurement (R).
P0  = blkdiag(R, 2^2 * eye(2));

% The initial estimate is the initial measurement and a random velocity
% error drawn from the initial covariance (just to make things
% interesting).
xh0 = [z0; x0(3:4) + covdraw(P0(3:4,3:4))]; % Initial estimate

% We'll linearize around a nominal state that's at terminal velocity.
vt = sqrt(9.81/cd);
xn = [0; 0; 0; -vt];

% Create the discrete state transition matrix from the continuous-time
% dynamics of the error state.
A  = [zeros(2), eye(2); zeros(2), -cd/m*vt * eye(2)]; % Continuous
F  = expm(A*dt);                                      % Discrete

% Create the observation matrix.
H  = [eye(2), zeros(2)];

% Create the constant input vector.
u  = f(xn, [0; 0]) - xn;
Fu = eye(4);
Hu = zeros(2, 4);

% Create the effective process noise for the filter (the effect of the
% process noise [acceleration] on the state update).
Fa = [0.5*dt^2 * eye(2); ... % Map acceleration to state
      dt * eye(2)];
Qe = Fa * Q * Fa.';          % Effective process noise matrix

% Create the true trajectory, measurements, and estimates.
t  = 0:dt:5;               % Time history
n  = length(t);            % Number of steps to take
x  = [x0, zeros(4, n-1)];  % True state history
q  = covdraw(Q, n-1);      % Process noise history (drawn up front)
z  = [z0, zeros(2, n-1)];  % Preallocation for the measurements
r  = covdraw(R, n);        % Measurement error history (drawn up front)
xh = [xh0, zeros(4, n-1)]; % Preallocation for the estimate history
P  = cat(3, P0, zeros(4, 4, n-1)); % Preallocation for the covariance

% Simulate the truth, create the noisy measurements, and run the filter for
% each step.
for k = 2:n
    x(:,k) = f(x(:,k-1), q(:,k-1));
    z(:,k) = x(1:2,k) + r(:,k);
    dx = xh(:,k-1) - xn; % Create the error state.
    [dx, P(:,:,k)] = lkf(dx, P(:,:,k-1), u, z(:,k), F, Fu, H, Hu, Qe, R);
    xh(:,k) = dx + xn;   % Rebuild the full state from the error state.
end

% Record the results with names for comparison to the EKF (run ekf_demo.m
% before this to compare with its results).
x_lkf  = x;
xh_lkf = xh;

%% Final plots

% Show the trajectory plot.
ell0 = ellipse(P0, xh0);
ell = ellipse(P(:,:,end), xh(:,end));
set(clf(figure(2)), 'Color', [1 1 1]);
ht = plot(x(1,:),    x(2,:), ...
          z(1,:),    z(2,:),  '.', ...
          xh(1,:),   xh(2,:), 'o', ...
          ell(1,:), ell(2,:));
set(ht(3:4), 'Color', 0.75 * [1 1 1]);
axis equal;
axis([x0(1)-1, x(1,end)+1, 0, x0(2)+1]);

%% Comparison
%
% Compare the results with the EKF (if its results are in the workspace).

if exist('xh_ekf', 'var')

    % States
    clf(figure(1));
    for k = 1:4
        subplot(4, 1, k);
        args = {t, x(k,:), ...
                t, xh_ekf(k,:), '--', ...
                t, xh_lkf(k,:), '--'};
        if k <= 2
            args = [args, {t, z(k,:), '.'}]; %#ok<AGROW>
        end
        plot(args{:});
        ylabel(sprintf('State %d', k));
    end
    xlabel('Time (s)');
    subplot(4, 1, 1);
    title('States over Time')
    legend('Truth', 'EKF', 'LKF', 'Measurements');

    % Errors
    clf(figure(3));
    for k = 1:4
        subplot(4, 1, k);
        args = {t, xh_ekf(k,:) - xh_lkf(k,:)};
        plot(args{:});
        ylabel(sprintf('State %d', k));
    end
    xlabel('Time (s)');
    subplot(4, 1, 1);
    title('Differences Between EKF and LKF');

end % if EKF results exist

%% Extended Kalman Filter Demonstration
%
% This script demonstrates using an extended Kalman filter to determine
% the position and velocity of a falling package for the article "How
% Kalman Filters Work".
%
% <http://www.anuncommonlab.com/articles/how-kalman-filters-work/>
%
% Copyright 2016 Tucker McClure @ An Uncommon Lab

%% Show the ambiguity of slope near the bounce.
%
% First, let's show that an EKF won't work on the ball problem (without
% some modifications) because the slope (and therefore the Jacobian) is
% undefined across the bounce.

% Run the ball simulation for a couple of seconds.
x0 = [0; 3; 1; 0];
[~, x, t] = propagate_ball(0, 2, x0);

% Find the index of the bounce.
bounce = find(x(2,:) == 0, 1);

% Create the dy/dt.
slope = diff(x(2,:)) ./ diff(x(1,:));

% Create lines corresponding to the slope before and after the bounce.
Dt  = 0.3;
dy1 = [-Dt * slope(bounce-1), 0];
dy2 = [0, Dt * slope(bounce)];

% Plot it.
clf(figure(2));
axis equal;
plot(x(1,1:2*bounce), x(2,1:2*bounce));
hold on;
plot(t(bounce) + [-Dt, 0], dy1);
plot(t(bounce) + [0, Dt],  dy2);
xlabel('Time [s]');
ylabel('Height [m]');
legend('Trajectory', 'Slope Before Bounce', 'Slope After Bounce', ...
       'Location', 'southeast');

%% System definition and simulation
%
% Instead of the bouncing ball, we'll use a parachuting package of coffee 
% filters instead. Let's define the system, set up the filter, and simulate
% the results over 10s.

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
Q = 0.5^2 * eye(2);  % Acceleration noise, 0.5m/s^2 std. dev.
R = 0.5^2 * eye(2);  % 0.5m standard deviation

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

% Propagation and observation functions for the filter.
ff = @(~, ~, x, ~) f(x, [0; 0]);
hf = @(~, x, ~) x(1:2);

% Create a function to calculate the Jacobian of the continuous dynamics 
% wrt the state.
dxdotdv = @(v) -cd/m*norm(v) * ((v*v.')/(v.'*v) + eye(2));
A       = @(x) [zeros(2), eye(2); ...
                zeros(2), dxdotdv(x(3:4))];

% Create a function to discretize the Jacobian of the continuous dynamics
% to obtain the Jacobian of the discrete propagation function.
F_fcn = @(~, ~, x, ~) expm(A(x) * dt);

% The observation matrix is constant, so create a function to return it.
H     = [eye(2), zeros(2)];
H_fcn = @(varargin) H;

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
    [xh(:,k), P(:,:,k)] = ekf(t(k-1), t(k), xh(:,k-1), P(:,:,k-1), ...
                              [], z(:,k), ff, hf, F_fcn, H_fcn, Qe, R);
end

% Record the results with names for comparison with the LKF (run lkf_demo.m
% after this one).
x_ekf  = x;
xh_ekf = xh;

%% Truth plot

% Make some axes limits to use repeatedly for consistency.
axes_limits = [x0(1)-1, x(1,end)+1, 0, x0(2)+1];

% Create the figure.
clf(figure(1));
plot(x(1,:), x(2,:));
title('Sample Trajectory of Dropped Package');
xlabel('x [m]');
ylabel('y [m]');
axis equal;
axis(axes_limits);

%% Propagation of state and covariance
% 
% Show the propagation of the first estimate in detail.

% We don't just want the updated state, we want lots of little points along
% the trajectory so you can see the curve of the propagation. We'll use 10
% little time steps to do this.
xht = [xh0, zeros(4, 10)];
for k = 2:11
    xht(:, k) = rk4step(fc, 0, xht(:, k-1), 0.1 * dt, [0; 0]);
end
xkkm1 = xht(:,end);

% Create the Jacobian and propagated covariance too.
Fkm1  = F_fcn(0, dt, xh0, []);
Pkkm1 = Fkm1 * P0 * Fkm1.' + Qe;

% Create the original and propagated ellipses.
ellkm1  = ellipse(P0, xh0);
ellkkm1 = ellipse(Pkkm1, xkkm1);

% Start over on the figure, zooming in right at the beginning.
clf(figure(1));

% Show the initial estimate.
hxh = plot(xh0(1), xh0(2), 'o', 'Color', 0.75 * [1 1 1]);
axis equal;
axis([xh0(1)-2 xh0(1)+2 xh0(2)-2 xh0(2)+2]);

% Add on the initial measurement and uncertainty.
hold on;
hz   = plot(z(1,2), z(2,2), '.');
hP0  = plot(ellkm1(1,:), ellkm1(2,:), 'Color', 0.75 * [1 1 1]);

% Label stuff.
legend([hxh, hP0, hz], ...
       'Initial Estimate', 'Initial Covariance', 'New Measurement');

%% Plot of propagated state

% Add on the ellipse for the propagated covariance.
set(hxh, 'XData', xkkm1(1), 'YData', xkkm1(2));
hxht = plot(xht(1,:), xht(2,:), 'Color', 0.75 * [1 1 1]);
legend([hxh, hP0, hz], ...
       'Propagated Estimate', 'Initial Covariance', 'New Measurement');

%% Plot of propated covariance

hPt = plot(ellkkm1(1,:), ellkkm1(2,:), ':', 'Color', 0.75 * [1 1 1]);
legend([hxh, hP0, hPt, hz], ...
       'Propagated Estimate', 'Initial Covariance', ...
       'Propagated Covariance', 'New Measurement');

%% Plot of Pxy and Pzz

% Delete the initial covariance ellipse.
delete(hP0);

% Calculate Pzz directly. Also, note that Pxy = P(1:2, 1:2) for our
% example.
Pzz = H * Pkkm1 * H.' + R;

% Add on the Pzz ellipse, and re-purpose the P ellipse for Pxy since it's
% the same thing.
ellPzz = ellipse(Pzz, xkkm1);
hPzz = plot(ellPzz(1,:), ellPzz(2,:), '--', 'Color', 0.75 * [1 1 1]);
legend([hxh, hPt, hPzz, hz], ...
       'Propagated Estimate', 'Propagated Covariance', ...
       'Innovation Covariance', 'New Measurement');

%% Corrected state and covariance

% We're done with the Pzz ellipse.
delete(hPzz);

% Create the corrected ellipse.
ellPk = ellipse(P(:,:,2), xh(:,2));

hxhkkm1 = plot(xkkm1(1), xkkm1(2), 'x', 'Color', 0.75 * [1 1 1]);
set(hxh, 'XData', xh(1,2), 'YData', xh(2,2));
set(hPt, 'XData', ellPk(1,:), 'YData', ellPk(2,:), 'LineStyle', '-');
legend([hxhkkm1, hxh, hPt, hz], ...
       'Predicted Estimate', 'Corrected Estimate', ...
       'Corrected Covariance', 'New Measurement');

%% Final plot and animation

% Draw the truth, measurements, and estimates as line plots.
clf(figure(1));
for k = 1:4
    subplot(4, 1, k);
    args = {t, x(k,:), ...
            t, xh(k,:), ':'};
    if k <= 2
        args = [args, {t, z(k,:), '.'}]; %#ok<AGROW>
    end
    plot(args{:});
end

% Show the trajectory plot, starting with the initial conditions.
ell0 = ellipse(P0, xh0);
set(clf(figure(2)), 'Color', [1 1 1]);
ht = plot(x(1,1),    x(2,1), ...
          z(1,1),    z(2,1),  '.', ...
          xh(1,1),   xh(2,1), 'o', ...
          ell0(1,:), ell0(2,:));
set(ht(3:4), 'Color', 0.75 * [1 1 1]); % Make estimate and covariance grey.
axis equal;
axis(axes_limits);

% Start the animated GIF.
make_gif = true;
if make_gif
    animation_name = fullfile('animations', 'ekf_demo_animation.gif');
    [A, map] = rgb2ind(frame2im(getframe()), 256);
    imwrite(A, map, animation_name, 'gif', ...
            'LoopCount', inf, ...
            'DelayTime', 2);
end

% Animate the descent.
tic();
for k = 1:length(t)
    
    % Set everything for this frame.
    ell = ellipse(P(:,:,k), xh(:,k));
    set(ht(1), 'XData', x(1,1:k),  'YData', x(2,1:k));
    set(ht(2), 'XData', z(1,1:k),  'YData', z(2,1:k));
    set(ht(3), 'XData', xh(1,1:k), 'YData', xh(2,1:k));
    set(ht(4), 'XData', ell(1,:),  'YData', ell(2,:));
    
    % Add to the animated GIF.
    if make_gif
        [A, map] = rgb2ind(frame2im(getframe()), 256);
        delay = dt;
        if k == length(t) % On the last frame, delay a while.
            delay = 2;
        end
        imwrite(A, map, animation_name, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime',  delay);
    end
    
    % Wait until it's time for the next frame.
    while toc() < t(k)
        pause(0.01);
    end
    
end % animation loop

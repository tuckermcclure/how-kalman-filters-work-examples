%% Sigma Point Filter Demonstration
%
% This script demonstrates using a sigma-point filter to determine the
% position and velocity of a bouncing ball for the article "How Kalman
% Filters Work".
%
% <http://www.anuncommonlab.com/articles/how-kalman-filters-work/>
%
% Copyright 2016 Tucker McClure @ An Uncommon Lab

%% Create the true system and show the initial filter state.

% Set the random number generator seed so the results are the same every
% time we run the script. (Comment out this line to see different results
% every time.)
rng(1);

% Initial true state, measurement noise covariance, and measurement
x0 = [0; 3; 1; 0];
R  = 0.5^2 * eye(2);
z0 = x0(1:2) + covdraw(R);

% Initial state estimate and covariance
xh0 = [z0; 1; 0];
P0  = blkdiag(R, 2^2 * eye(2));

% Calculate the whole true trajectory.
[~, x, t] = propagate_ball(0, 10, x0);

% Prepare the figure.
set(clf(figure(1)), 'Color', [1 1 1]);
axis equal;
axis([-1 11 0 5]);
xlabel('x [m]');
ylabel('y [m]');
hold on;

% Draw the 3-sigma boundary for the uncertainty.
ell = ellipse(P0, xh0);
hP  = plot(ell(1,:), ell(2,:), 'Color', 0.75 * [1 1 1]);

% Add particles for comparison only.
X  = bsxfun(@plus, covdraw(P0, 1000), xh0);
hX = plot(X(1,:), X(2,:), '.', 'Color', 0.75 * [1 1 1]);

% Add on the initial measurement.
hz = plot(z0(1), z0(2), 'o', 'Color', [0.85 0.325 0.098]);

% Start the legend.
legend([hz hP hX], 'Measurement', '3\sigma boundary', 'Particles');

%% Add the inital sigma points.

% Get rid of the particles in favor of the sigma points.
delete(hX); hX = [];

% A function to create the (visualization-only) sigma points from the 
% covariance and state. These are different from the real sigma points in
% that they don't include process or measurement noise and use a very large
% scaling factor (1), whereas the real filter will have the sigma points
% spaced too closely to see clearly.
sigma_points = @(P, x) ...
    bsxfun(@plus, ...
           1 * chol(P, 'lower') * [[0; 0; 0; 0], eye(4), -eye(4)], ...
           x);

% Initial sigma points
X  = sigma_points(P0, xh0);
nX = size(X, 2);

% Add the initial sigma points.
hS = plot(X(1,:), X(2,:), '.', 'Color', 0.75 * [1 1 1]);

% Update the legend.
legend([hz hP hS], 'Measurement', '3\sigma boundary', 'Sigma Points');

%% Propagate the truth and draw the measurement.

% Manually propagate the truth and take a measurement.
tk = 0.5;
xk = propagate_ball(0, tk, x0);
zk = xk(1:2) + covdraw(R);
set(hz, 'XData', [get(hz, 'XData'), zk(1)], ...
        'YData', [get(hz, 'YData'), zk(2)]);

%% Propagate the sigma points, drawing their trajectories.

% Get trajectories for each ball.
xt = cell(1, nX);
ht = zeros(1, nX);
for k = 1:nX
    [~, xt{k}] = propagate_ball(0, tk, X(:,k));
    ht(k) = plot(xt{k}(1,:), xt{k}(2,:), 'Color', 0.75 * [1 1 1]);
end

%% Run the UKF and show the predicted stuff.

% Set up the propagation and measurement functions for the UKF.
f = @propagate_ball;
h = @(t, x, u, r) x(1:2) + r;

% Our problem doesn't have process noise.
Q = [];

% Run the filter.
[xhk, Pk, xkkm1, Pkkm1, Pxy, Pyy, K] = ukf(0, tk, xh0, P0, [], zk, ...
                                           f, h, Q, R); % , 0.01, 2, -1

% Show the predicted state and propagated covariance.
hxh = scatter(xkkm1(1), xkkm1(2), 100, 0.75*[1 1 1], 'o', 'filled');
ell = ellipse(Pkkm1, xhk);
set(hP, 'XData', ell(1,:), 'YData', ell(2,:));

%% Add on P_xy and P_yy.

% Hide the estimate covariance for a moment.
set(hP, 'Visible', 'off');

% Add on the state-innovation covariance and the innovation covariance.
ell = ellipse(Pxy, xhk);
hPxy = plot(ell(1,:), ell(2,:), ':', 'Color', 0.75*[1 1 1]);
ell = ellipse(Pyy, xhk);
hPyy = plot(ell(1,:), ell(2,:), '--', 'Color', 0.75*[1 1 1]);
legend([hz hS hPxy hPyy], ...
       'Measurement', 'Sigma Points', ...
       '3\sigma of P_{xy}', '3\sigma of P_{yy}');

%% Show the results.

% Get rid of the intermediate stuff.
delete(hPxy);
delete(hPyy);
set(hP, 'Visible', 'on');

% Show the corrected state estimate and covariance.
set(hxh, 'XData', xhk(1), 'YData', xhk(2));
ell = ellipse(Pk, xhk);
set(hP, 'XData', ell(1,:), 'YData', ell(2,:));

% Update the legend.
legend([hz hP hS hxh], ...
       'Measurement', '3\sigma boundary', 'Sigma Points', 'Estimate');

%% Run out the filter for 10s.

% Remove the old sigma points.
delete(hS); hS = [];

% Add on the truth.
hx = plot(xk(1), xk(2), 'x', 'Color', [0 0.4470 0.7410]);

% Update the legend.
legend([hz hP hxh hx], ...
       'Measurement', '3\sigma boundary', 'Estimate', 'Truth');

% Start the animated GIF.
make_gif = true;
if make_gif
    animation_name = fullfile('animations', 'sigma_point_demo_animation.gif');
    [A, map] = rgb2ind(frame2im(getframe()), 256);
    imwrite(A, map, animation_name, 'gif', ...
            'LoopCount', inf, ...
            'DelayTime', 2);
end

% Set the time step for the simulation.
dt = 0.1;

% Loop until we reach 10s, animating the results.
tic();
t0 = tk;
for tk = tk+dt:dt:10
    
    % Update truth and create noisy measurement.
    xk = propagate_ball(tk-dt, tk, xk);
    zk = xk(1:2) + covdraw(R);
    
    % Create little trajectories for sigma points (just for the plot).
    X = sigma_points(Pk, xhk);
    for k = 1:size(X,2)
        [~, xt{k}] = propagate_ball(tk-dt, tk, X(:,k));
        set(ht(k), 'XData', xt{k}(1,:), ...
                   'YData', xt{k}(2,:));
    end

    % Run the filter.
    [xhk, Pk] = ukf(tk-dt, tk, xhk, Pk, [], zk, f, h, Q, R); % , 0.01, 2, -1
    
    % Draw everything.
    ell = ellipse(Pk, xhk);
    set(hP,  'XData', ell(1,:), ...
             'YData', ell(2,:));
    set(hz,  'XData', [get(hz, 'XData') zk(1)], ...
             'YData', [get(hz, 'YData') zk(2)]);
    set(hxh, 'XData', [get(hxh, 'XData') xhk(1)], ...
             'YData', [get(hxh, 'YData') xhk(2)]);
    set(hx,  'XData', [get(hx, 'XData') xk(1)], ...
             'YData', [get(hx, 'YData') xk(2)]);
        
    % Add to the animated GIF.
    if make_gif
        [A, map] = rgb2ind(frame2im(getframe()), 256);
        delay = dt;
        if tk == 10
            delay = 2;
        end
        imwrite(A, map, animation_name, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', delay);
    end
    
    % Wait until it's time for the next frame.
    while toc() + t0 < tk + dt
        pause(0.01);
    end

end % simulation loop

% Plop full truth on top of it.
htt = plot(x(1,:), x(2,:), 'Color', [0 0.4470 0.7410]);

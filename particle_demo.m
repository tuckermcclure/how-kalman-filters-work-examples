%% Particle Filter Demonstration
%
% This script demonstrates using a particle filter to determine the
% position and velocity of a bouncing ball for the article "How Kalman
% Filters Work".
%
% <http://www.anuncommonlab.com/articles/how-kalman-filters-work/>
%
% Copyright 2016 Tucker McClure @ An Uncommon Lab

%% Set up and start the figure.

% Set the random number generator seed so the results are the same every
% time we run the script. (Comment out this line to see different results
% every time.)
rng(1);

% Initial true state, measurement noise covariance, and measurement
x0 = [0; 3; 1; 0];
R  = 0.25^2 * eye(2);
z0 = x0(1:2) + covdraw(R);

% Initial state estimate and covariance. The initial particles will be
% centered on this initial estimate and distributed according to the
% covariance. The initial particle weights are all equal (1/N).
xh0 = [z0; 1; 0];
P0  = blkdiag(R, 2^2 * eye(2));
nX  = 100;
X   = bsxfun(@plus, covdraw(P0, nX), xh0);
w   = 1/nX * ones(1, nX);

% Calculate the whole true trajectory.
[~, x, t] = propagate_ball(0, 10, x0);

% Create a function to calculate the colors for the particles according to
% their weights. This looks pretty good.
colors = @(w) 1 - max(0.95*w/max(w).^0.85, 0.05);

% Prepare the figure.
set(clf(figure(1)), 'Color', [1 1 1]);
colormap('gray');
axis equal;
axis([-1 11 0 5]);
caxis([0 1]);
xlabel('x [m]');
ylabel('y [m]');
hold on;

% Add the initial particles.
hX = scatter(X(1,:), X(2,:), [], colors(w), '.');

% Measurement
hz = plot(z0(1), z0(2), 'o', 'Color', [0.85 0.325 0.098]);

% Start the legend.
legend([hz hX], 'Measurement', 'Particles');
% colorbar();

%% Propagate the truth and draw the measurement.

% Manually propagate the truth and take a noisy measurement.
tk = 1;
xk = propagate_ball(0, tk, x0);
zk = xk(1:2) + covdraw(R);
set(hz, 'XData', [get(hz, 'XData'), zk(1)], ...
        'YData', [get(hz, 'YData'), zk(2)]);

%% Propagate the particles, drawing their trajectories.

% Get trajectories for each ball.
xt = cell(1, nX);
Xn = zeros(size(X));
ht = zeros(1, nX);
line_colors = colors(w).' * [1 1 1];
for k = 1:nX
    [Xn(:,k), xt{k}] = propagate_ball(0, 1, X(:,k));
    ht(k) = plot(xt{k}(1,:), xt{k}(2,:), 'Color', line_colors(k,:));
end
set(hX, 'XData', Xn(1,:), 'YData', Xn(2,:), 'CData', colors(w));

%% Show distance to each particle.

hd = zeros(1, nX);
for k = 1:nX
    hd(k) = plot([xt{k}(1,end) zk(1)], [xt{k}(2,end), zk(2)], ...
                 'Color', [0.85 0.325 0.098]);
end

%% Use the bootstrap filter to update the weights.

% Get rid of the error lines.
delete(hd);

% Create the propagation and observation functions.
f = @propagate_ball;
h = @(t, x, u) x(1:2);

% Create a function to make a random process noise draw. Since we don't use
% process noise for the bouncing ball, this is pretty easy:
d = @(varargin) [];

% Create a function to determine how likely a measurement error is. Since
% our measurement errors are Gaussian, we'll use the probability density
% function for a multivariate Gaussian distribution. (Note, since the
% "probabilities" of each particle are only important relative to each
% other, we can drop the constant scaling factor of the PDF.)
invR = inv(R);
p = @(t, dz, varargin) exp(-0.5 * dz.' * invR * dz); %#ok<MINV>

% Run the bootstrap filter for the first time step.
[xh, X, w, wt] = bf(0, tk, ...  % Last and current times
                    X, w, ...   % Last particles and weights
                    [], ...     % Input vector (we don't use one)
                    zk, ...     % Current measurement
                    f, h, ...   % Propagation and measurement functions
                    d, ...      % Function to draw random process noise
                    p, ...      % Fcn to determine prob. of meas. error
                    [], ...     % Tuning parameter (let it use a default)
                    false);     % Don't resample and regularize (see below)

% Update the particle plot.
set(hX, 'XData', X(1,:), 'YData', X(2,:), 'CData', colors(wt));

% Update the trajectories.
line_colors = colors(wt).' * [1 1 1];
for k = 1:nX
    set(ht(k), 'Color', line_colors(k,:));
end

%% Show weighted average.

hxh = scatter(xh(1), xh(2), 100, [0.25 0.75 0.25], 'o', 'filled');
legend([hz hX hxh], 'Measurement', 'Particles', 'Estimated State');

%% Propagate again.

% Propagate the truth and create a noisy measurement.
tkm1 = tk;
tk   = 2;
xk = propagate_ball(tkm1, tk, xk);
zk = xk(1:2) + covdraw(R);

% Create the little trajectories for each particle.
for k = 1:nX
    [Xn(:,k), xt{k}] = propagate_ball(tkm1, tk, X(:,k));
end

% Update the measurement.
set(hz, 'XData', [get(hz, 'XData'), zk(1)], ...
        'YData', [get(hz, 'YData'), zk(2)]);

% Update the particles.
set(hX, 'XData', Xn(1,:), 'YData', Xn(2,:), 'CData', colors(wt));

% Update their trajectories.
line_colors = colors(wt).' * [1 1 1];
for k = 1:nX
    set(ht(k), 'XData', xt{k}(1,:), 'YData', xt{k}(2,:), ...
        'Color', line_colors(k,:));
end

%% Update weights again.

% Just get the updated weights. We don't want the new particles, because
% we're going to show the resampled and regularized particles next.
[~, ~, ~, wt] = bf(tkm1, tk, X, w, [], zk, f, h, d, p, [], false);
% draw_particles(Xn, wt);

% Update the particles.
set(hX, 'CData', colors(wt));

% Update their trajectories.
line_colors = colors(wt).' * [1 1 1];
for k = 1:nX
    set(ht(k), 'Color', line_colors(k,:));
end

%% Resample and regularize.

% Run the filter again, this time resampling and regularizing the
% particles.
[xh, X, w, wt] = bf(tkm1, tk, X, w, [], zk, f, h, d, p, [], true);

% Show the new estimate.
set(hxh, 'XData', [get(hxh, 'XData'), xh(1)], ...
         'YData', [get(hxh, 'YData'), xh(2)]);

% Update the particles.
set(hX, 'XData', X(1,:), 'YData', X(2,:), 'CData', colors(w));

% Hide the trajectories; they no longer connect to the particles.
set(ht, 'Visible', 'off');

% return;

%% Run out the filter for 10s.

% Turn the trajectories back on.
set(ht, 'Visible', 'on');

% Truth
hx = plot(xk(1), xk(2), 'x', 'Color', [0 0.4470 0.7410]);
legend([hz hX hxh hx], 'Measurements', 'Particles', 'Estimated', 'Truth');

% Start the animated GIF.
make_gif = true;
if make_gif
    animation_name = fullfile('animations', 'particle_demo_animation.gif');
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
    
    % Create little trajectories for particles (just for the plot).
    for k = 1:nX
        [~, xt{k}] = propagate_ball(tk-dt, tk, X(:,k));
    end

    % Resample on every other time step.
    resample = mod(tk, 0.2) < eps;
    
    % Run the filter.
    [xhk, X, w, wt] = bf(tk-dt, tk, X, w, [], zk, f, h, d, p, [], ...
                         resample);

    % Draw everything.
    set(hX,  'XData', X(1,:), 'YData', X(2,:), 'CData', colors(w));
    line_colors = colors(wt).' * [1 1 1];
    for k = 1:nX
        set(ht(k), 'XData', xt{k}(1,:), ...
                   'YData', xt{k}(2,:), ...
                   'Color', line_colors(k,:));
    end
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
htt = plot(x(1, :), x(2, :), 'Color', [0 0.4470 0.7410]);

%% Example multimodal probability distribution

% This simple function will create the matrix square root of a random,
% valid covariance matrix.
randsqrtcov = @()   orth(randn(2)) ...      % Random rotation matrix
                  * diag(sqrt(rand(1, 2))); % Random eigenvalues

% Draw 3 random Gaussian distributions about 3 distinct means.
rng(2);
n = 500;
figure(2);
a = [randsqrtcov() * randn(2, n) + repmat(3 * randn(2,1), 1, n), ...
     randsqrtcov() * randn(2, n) + repmat(3 * randn(2,1), 1, n), ...
     randsqrtcov() * randn(2, n) + repmat(3 * randn(2,1), 1, n)];
plot(a(1,:), a(2,:), '.', 'Color', 0.75*[1 1 1]);
title('Example of Particles for a Multimodal Probability Distribution');
xlabel('State 1');
ylabel('State 2');

%% Set up and start the figure.

% Setup
clear all; %#ok<CLALL>
chill = @(x) drawnow();
% chill = @(x) pause(x);
rng(1);

% Initial system
x0 = [0; 3; 1; 0];
P0 = diag([0.25 0.25 3 3]);
R  = 0.25 * eye(2);

% Initial measurement and state estimate
z0  = x0(1:2) + mnddraw(R);
xh0 = [z0; 1; 0];

% Calculate the whole true trajectory.
[~, x] = propagate(0, 10, x0);

% Prepare the figure.
set(clf(figure(1)), 'Color', [1 1 1]);
axis equal;
axis([-1 11 0 5]);
xlabel('x [m]');
ylabel('y [m]');
hold on;

% Tools to draw ellipses
circle  = feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));
ellipse = @(P,x) bsxfun(@plus, 3 * chol(P(1:2,1:2), 'lower') * circle, x(1:2));

% 3-sigma boundary
ell = ellipse(P0, xh0);
hP  = plot(ell(1,:), ell(2,:), 'Color', 0.75 * [1 1 1]);

% Particles
X  = bsxfun(@plus, mnddraw(P0, 1000), xh0);
hX = plot(X(1,:), X(2,:), '.', 'Color', 0.75 * [1 1 1]);

% Measurement
hz = plot(z0(1), z0(2), 'o', 'Color', [0.85 0.325 0.098]);

% Start the legend.
legend([hz hP hX], 'Measurement', '3\sigma boundary', 'Particles');

%% Add the inital sigma points.

% Get rid of the particles in favor of the sigma points.
delete(hX); hX = [];

% Sigma point scaling factor
a = 1;

% A function to create the sigma points from the covariance and state
sigma_points = @(P, x) ...
    bsxfun(@plus, ...
           a * chol(P, 'lower') * [[0; 0; 0; 0], eye(4), -eye(4)], ...
           x);

% Initial sigma points
S  = sigma_points(P0, xh0);
nS = size(S, 2);

% Add the initial sigma points.
hS = plot(S(1,:), S(2,:), '.', 'Color', 0.75 * [1 1 1]);

% Update the legend.
legend([hz hP hS], 'Measurement', '3\sigma boundary', 'Sigma Points');

chill(0.5);

%% Propagate the truth and draw the measurement.

% Manually propagate the truth and take a measurement.
tk = 0.5;
xk = propagate(0, tk, x0);
zk = xk(1:2) + mnddraw(R);
plot(zk(1), zk(2), 'ro');

chill(0.5);

%% Propagate the sigma points, drawing their trajectories.

% Get trajectories for each ball.
xt = cell(1, nS);
ht = zeros(1, nS);
for k = 1:nS
    [~, xt{k}] = propagate(0, tk, S(:,k));
    ht(k) = plot(xt{k}(1,:), xt{k}(2,:), 'Color', 0.75 * [1 1 1]);
end

chill(0.5);

%% Run the UKF and show the predicted stuff.

f = @propagate;
h = @(~, x, varargin) x(1:2);
Q = zeros(4);
[xhk, Pk, xkkm1, Pkkm1, Pxy, Pyy, K] = ukfan(0, tk, xh0, P0, [], zk, f, h, Q, R, 0.01, 2, -1);

% f = @(t0,tf,x,u,q) propagate(t0,tf,x) + q;
% h = @(~, x, ~, r) x(1:2) + r;
% [xhk, Pk] = ukf(0, tk, xh0, P0, [], zk, f, h, Q, R, 0.01, 2, -1);

hxh = scatter(xkkm1(1), xkkm1(2), 100, 0.75*[1 1 1], 'o', 'filled');
ell = ellipse(Pkkm1, xhk);
set(hP, 'XData', ell(1,:), 'YData', ell(2,:));

%% Add on P_xy and P_yy.

set(hP, 'Visible', 'off');

ell = ellipse(Pxy, xhk);
hPxy = plot(ell(1,:), ell(2,:), ':', 'Color', 0.75*[1 1 1]);
ell = ellipse(Pyy, xhk);
hPyy = plot(ell(1,:), ell(2,:), '--', 'Color', 0.75*[1 1 1]);
legend([hz hS hPxy hPyy], 'Measurement', 'Sigma Points', '3\sigma of P_{xy}', '3\sigma of P_{yy}');

%% Show the results.

delete(hPxy);
delete(hPyy);
set(hP, 'Visible', 'on');

% draw_particles(sigma_points(Pk, xh), w);
set(hxh, 'XData', xhk(1), 'YData', xhk(2));
ell = ellipse(Pk, xhk);
set(hP, 'XData', ell(1,:), 'YData', ell(2,:));

% Update the legend.
legend([hz hP hS hxh], 'Measurement', '3\sigma boundary', 'Sigma Points', 'Estimate');

%% Run out the filter for 10s.

% Remove the old sigma points.
delete(hS); hS = [];

% Add on the truth.
hx = plot(xk(1), xk(2), 'x', 'Color', [0    0.4470    0.7410]);%[0.4660 0.6740 0.1880]); %[0.4940 0.1840 0.5560]);

% Update the legend.
legend([hz hP hxh hx], 'Measurement', '3\sigma boundary', 'Estimate', 'Truth');

% Start the animated GIF.
animation_name = fullfile('..', 'jade', 'img', 'sigma_point_demo_animation.gif');
[A, map] = rgb2ind(frame2im(getframe()), 256);
imwrite(A, map, animation_name, 'gif', 'LoopCount', inf, 'DelayTime', 2);

dt = 0.1;
for tk = tk+dt:dt:10
    
    % Update truth and measure.
    xk = propagate(tk-dt, tk, xk);
    zk = xk(1:2) + mnddraw(R);
    
    % Create trajectories for sigma points.
    S = sigma_points(Pk, xhk);
    for k = 1:size(S,2)
        [~, xt{k}] = propagate(tk-dt, tk, S(:,k));
        set(ht(k), 'XData', xt{k}(1,:), ...
                   'YData', xt{k}(2,:));
    end

    % Run the filter.
    [xhk, Pk] = ukfan(tk-dt, tk, xhk, Pk, [], zk, f, h, Q, R, 0.01, 2, -1);
    
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
    
    chill(dt);
    
    % Add to the animated GIF.
    [A, map] = rgb2ind(frame2im(getframe()), 256);
    delay = dt;
    if tk == 10
        delay = 2;
    end
    imwrite(A, map, animation_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay);

end

% Plop truth on top of it.
htt = plot(x(1,:), x(2,:), 'Color', [0    0.4470    0.7410]);%[0.4660 0.6740 0.1880]); %[0.4940 0.1840 0.5560]);

%%

% Setup
clear all; %#ok<CLALL>
rng(1);

% Calculate the whole true trajectory.
x0 = [0; 3; 1; 0];
[~, x] = propagate(0, 2, x0);

bounce = find(x(2,:) == 0, 1);
slope = diff(x(2,:)) ./ diff(x(1,:));
dx = 0.3;
dx1 = [x(1:2,bounce) - dx*[1; slope(bounce-1)], x(1:2,bounce) + 0*dx*[1; slope(bounce-1)]];
dx2 = [x(1:2,bounce) - 0*dx*[1; slope(bounce)],   x(1:2,bounce) + dx*[1; slope(bounce)]];

% Prepare the figure.
set(clf(figure(1)), 'Color', [1 1 1]);
axis equal;
plot(x(1,:), x(2,:));
hold on;
plot(dx1(1,:), dx1(2,:));
plot(dx2(1,:), dx2(2,:));
% axis([-1 3 0 5]);
xlabel('x [m]');
ylabel('y [m]');
legend('Trajectory', 'Slope Before Bounce', 'Slope After Bounce');

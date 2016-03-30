%% Set up and start the figure.

% Setup
clear all; %#ok<CLALL>
rng(1);

% Calculate the whole true trajectory.
x0 = [0; 3; 1; 0];
[~, x, t] = propagate(0, 2, x0);

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

%%

rng(1);

% Start at a bounce.
x0 = [0; 0; 2; 10];
P0 = diag([0.25 0.25 3 3]);
R  = 0.25 * eye(2);

% Initial measurement and state estimate
z0  = x0(1:2) + mnddraw(R);
xh0 = [z0; 2; 11];

% Calculate the whole true trajectory.
[~, x, t] = propagate(0, 2, x0);

% Propagate the initial estimate.
[~, xh, th] = propagate(0, 0.1, xh0);

% Propagate the state and covariance.

% Prepare the figure.
set(clf(figure(1)), 'Color', [1 1 1]);
axis equal;
axis([-1 5 0 6]);
xlabel('x [m]');
ylabel('y [m]');
hold on;

% Tools to draw ellipses
circle  = feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));
ellipse = @(P,x) bsxfun(@plus, 3 * chol(P(1:2,1:2), 'lower') * circle, x(1:2));

% 3-sigma boundary
ell = ellipse(P0, xh0);
hP  = plot(ell(1,:), ell(2,:), 'Color', 0.75 * [1 1 1]);

hx = plot(x(1,:), x(2,:));
hxh = plot(xh(1,:), xh(2,:), 'Color', 0.75 * [1 1 1]);

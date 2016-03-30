%%
chill = @(x) drawnow();
% chill = @(x) pause(x);

%% Set up the initial particles.
rng(1);
% capture('setup', 'pf');

x0 = [0; 3; 1; 0];
nX = 100;
P  = diag([0.25^2 0.25^2 3 3]);
R  = 0.25^2 * eye(2);

z0 = x0(1:2);
X  = bsxfun(@plus, mnddraw(P, nX), [z0; 0; 0]);
w  = 1/nX * ones(1, nX);

draw_particles('reset');

%% Calculate the whole true trajectory.

[~, x] = propagate(0, 10, x0);

%% Draw the inital particles and measurement.

% Prepare the figure.
figure(1);
clf();
colormap('gray');
axis([-1 11 0 5]);
xlabel('x [m]');
ylabel('y [m]');
hold on;
axis equal;

% Add the initial particles.
draw_particles(X, w);

% Add on the measurement.
h = plot(z0(1), z0(2), 'ro');

chill(0.5);
% capture();

%% Propagate the truth and draw the measurement.

% Manually propagate the truth and take a measurement.
xk = propagate(0, 1, x0);
zk = xk(1:2) + mnddraw(R);
plot(zk(1), zk(2), 'ro');

chill(0.5);

%% Propagate the particles, drawing their trajectories.

% Get trajectories for each ball.
xt = cell(1, nX);
Xn = zeros(size(X));
for k = 1:nX
    [Xn(:,k), xt{k}] = propagate(0, 1, X(:,k));
end
draw_particles(Xn, w, xt);

chill(0.5);

%% Show distance to each particle.

hd = zeros(1, nX);
for k = 1:nX
    hd(k) = plot([xt{k}(1, end) zk(1)], [xt{k}(2, end), zk(2)], 'r');
end

%% Use the bootstrap filter to update the weights.

delete(hd);

invR = inv(R);
f = @propagate;
y = @(t, x, varargin) x(1:2);
d = @(varargin) [];
p = @(t, dz, varargin) exp(-0.5 * dz.' * invR * dz); %#ok<MINV>
[xh, X, w, wt] = bf(0, 1, X, w, [], zk, f, y, d, p, [], false);
draw_particles(X, wt);

chill(0.5);

%% Show weighted average.

scatter(xh(1), xh(2), 100, [0.05 0.95 0.05], 'o', 'filled');

chill(0.5);

%% Propagate again.

xk = propagate(1, 2, xk);
zk = xk(1:2) + mnddraw(R);

for k = 1:nX
    [Xn(:,k), xt{k}] = propagate(0, 1, X(:,k));
end
draw_particles(Xn, w, xt);

plot(zk(1), zk(2), 'ro');

chill(0.5);

%% Update weights again.

[~, ~, ~, wt] = bf(1, 2, X, w, [], zk, f, y, d, p, [], false);
draw_particles(Xn, wt);

chill(0.5);

%% Resample and regularize.

[xh, X, w, wt] = bf(1, 2, X, w, [], zk, f, y, d, p, [], true);
draw_particles(X, w, []);
scatter(xh(1), xh(2), 100, [0.05 0.95 0.05], 'o', 'filled');

chill(0.5);

%% Run out the filter for 10s.

% Start the animated GIF.
animation_name = fullfile('jade', 'img', 'particle_demo_animation.gif');
[A, map] = rgb2ind(frame2im(getframe()), 256);
imwrite(A, map, animation_name, 'gif', 'LoopCount', inf, 'DelayTime', 2);

dt = 0.1;
axis([-1 11 0 4]);
for t = 2+dt:dt:10
    
    xk = propagate(t-dt, t, xk);
    zk = xk(1:2) + mnddraw(R);
    
    for k = 1:nX
        [~, xt{k}] = propagate(t-dt, t, X(:,k));
    end

    [xh, X, w, wt] = bf(t-dt, t, X, w, [], zk, f, y, d, p, [], true);

    % Draw everything.
    hp = draw_particles([], wt, xt); % Particles (w/o velocity arrows)
    hm = plot(zk(1), zk(2), 'ro');  % Measured
    he = scatter(xh(1), xh(2), 100, [0.05 0.95 0.05], 'o', 'filled'); % Estimated
    ht = plot(xk(1), xk(2), 'bo');  % True
    legend([hp, hm, he, ht], 'Particles', 'Measured', 'Estimated', 'Truth');
    
    chill(0.15);
    
    % Add to the animated GIF.
    [A, map] = rgb2ind(frame2im(getframe()), 256);
    delay = dt;
    if t == 10
        delay = 2;
    end
    imwrite(A, map, animation_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay);

end

% Plop truth on top of it.
ht = plot(x(1, :), x(2, :), 'b');

%% Multimodal probability distribution

rng(9);
n = 100;
figure(2);
a = [randcov(2, [], [], 'sqrt') * randn(2, n) + repmat(3 * randn(2,1), 1, n), ...
     randcov(2, [], [], 'sqrt') * randn(2, n) + repmat(3 * randn(2,1), 1, n), ...
     randcov(2, [], [], 'sqrt') * randn(2, n) + repmat(3 * randn(2,1), 1, n)];
scatter(a(1,:), a(2,:));
title('Example of Particles for a Multi-Modal Probability Distribution');
xlabel('State 1');
ylabel('State 2');

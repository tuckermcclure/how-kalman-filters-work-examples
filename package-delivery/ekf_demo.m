%% Set up and start the figure.

addpath('../bouncing-ball');

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

rmpath('../bouncing-ball');

%% System definition and simulation

% cc;
clc;
rng(8);

% System
dt = 0.1;
cd = 4.4;
m  = 1;
Fa = [0.5*dt^2 * eye(2); dt * eye(2)];
Q  = 0.5^2 * eye(2);
Qe = Fa * Q * Fa.';
R  = 0.25^2 * eye(2);
g  = [0; -9.81];
fc = @(t, x, q) [x(3:4); g - cd/m * norm(x(3:4)) * x(3:4) + q];
f  = @(x, q) rk4step(fc, 0, x, dt, q);

% Initial conditions
x0  = [-1; 8; 1; 0];
P0  = bdiag(R, 1^2 * eye(2));
z0  = x0(1:2) + mnddraw(R);
xh0 = [z0; x0(3:4) + mnddraw(P0(3:4,3:4))];

% Create the filter options.
ff = @(~, ~, x, ~) f(x, [0; 0]);
hf = @(~, x, ~) x(1:2);
A  = @(x) [zeros(2), eye(2); zeros(2), -cd/m*norm(x(3:4)) * ((x(3:4)*x(3:4).')/(x(3:4).'*x(3:4)) + eye(2))];
F  = @(~, ~, x, ~) expm(A(x) * dt);
H  = [eye(2), zeros(2)];
options = kffoptions('f', ff, ...
                     'F_km1_fcn', F, ...
                     'Q_km1_fcn', Qe, ...
                     'h', hf, ...
                     'H_k_fcn', H, ...
                     'R_k_fcn', R);

% Create the true trajectory, measurements, and estimates.
t  = 0:dt:5;
n  = length(t);
x  = [x0, zeros(4, n-1)];
q  = mnddraw(Q, n-1);
z  = [z0, zeros(2, n-1)];
r  = mnddraw(R, n);
xh = [xh0, zeros(4, n-1)];
P  = cat(3, P0, zeros(4, 4, n-1));
for k = 2:n
    x(:,k) = f(x(:,k-1), q(:,k-1));
    z(:,k) = x(1:2,k) + r(:,k);
    [xh(:,k), P(:,:,k)] = kff(xh(:,k-1), P(:,:,k-1), z(:,k), options);
end

% Tools to draw ellipses
circle  = feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));
ellipse = @(P,x) bsxfun(@plus, 3 * chol(P(1:2,1:2), 'lower') * circle, x(1:2));

% Record the results with names for comparison with the LKF.
x_ekf = x; xh_ekf = xh;

%% Truth

clf(figure(1));
plot(x(1,:), x(2,:));
title('Sample Trajectory of Dropped Package');
xlabel('x [m]');
ylabel('y [m]');
axis equal;
axis([x0(1)-1, x(1,end)+1, 0, x0(2)+1]);

%% Propagation of state and covariance

xht = [xh0, zeros(4, 10)];
for k = 2:11
    xht(:, k) = rk4step(fc, 0, xht(:, k-1), 0.1 * dt, [0; 0]);
end
xkkm1 = xht(:,end);
Fkm1  = F(0, dt, x0, []);
Pkkm1 = Fkm1 * P0 * Fkm1.' + Qe;

ellkm1  = ellipse(P0, xh0);
ellkkm1 = ellipse(Pkkm1, xkkm1);

clf(figure(1));
hxh  = plot(xh0(1), xh0(2), 'o', 'Color', 0.75 * [1 1 1]);
axis equal;
axis([xh0(1)-1 xh0(1)+1 xh0(2)-1 xh0(2)+1]);
hold on;
hz   = plot(z(1,2), z(2,2), '.');
hP0  = plot(ellkm1(1,:), ellkm1(2,:), 'Color', 0.75 * [1 1 1]);
legend([hxh, hP0, hz], 'Initial Estimate', 'Initial Covariance', 'New Measurement');

%%

set(hxh, 'XData', xkkm1(1), 'YData', xkkm1(2));
hxht = plot(xht(1,:), xht(2,:), 'Color', 0.75 * [1 1 1]);
legend([hxh, hP0, hz], 'Propagated Estimate', 'Initial Covariance', 'New Measurement');

%%
hPt  = plot(ellkkm1(1,:), ellkkm1(2,:), ':', 'Color', 0.75 * [1 1 1]);
legend([hxh, hP0, hPt, hz], 'Propagated Estimate', 'Initial Covariance', 'Propagated Covariance', 'New Measurement');

%% Pxy and Pzz

delete(hP0);

Pzz = H * Pkkm1 * H.' + R;
ellPzz = ellipse(Pzz, xkkm1);
hPzz = plot(ellPzz(1,:), ellPzz(2,:), '--', 'Color', 0.75 * [1 1 1]);
legend([hxh, hPt, hPzz, hz], 'Propagated Estimate', 'Propagated Covariance', 'Innovation Covariance', 'New Measurement');

%% Corrected state and covariance

delete(hPzz);

K = Pkkm1 * H.' / Pzz;
xk = xkkm1 + K * (z(:,2) - xkkm1(1:2));
Pk = Pkkm1 - K * H * Pkkm1;
ellPk = ellipse(Pk, xk);

set(hxh, 'XData', xk(1), 'YData', xk(2));
set(hPt, 'XData', ellPk(1,:), 'YData', ellPk(2,:), 'LineStyle', '-');
legend([hxh, hPt, hz], 'Corrected Estimate', 'Corrected Covariance', 'New Measurement');

%% Final plot and animation

% Draw the truth, measurements, and estimates as line plots.
clf(figure(1));
% hl = zeros(1, 4);
for k = 1:4
    subplot(4, 1, k);
    args = {t, x(k,:), ...
            t, xh(k,:), ':'};
    if k <= 2
        args = [args, {t, z(k,:), '.'}]; %#ok<AGROW>
    end
    plot(args{:});
end
% rowsubplot(@plot, t, x, t, xh, t, z, {{}, {':'}, {'.'}});
% rowsubplot(@plot, {t, x}, {t, xh, ':'}, {t, z, '.'});

% Show the trajectory plot.
ell0 = ellipse(P0, xh0);
set(clf(figure(2)), 'Color', [1 1 1]);
ht = plot(x(1,1),    x(2,1), ...
          z(1,1),    z(2,1),  '.', ...
          xh(1,1),   xh(2,1), 'o', ...
          ell0(1,:), ell0(2,:));
set(ht(3:4), 'Color', 0.75 * [1 1 1]);
axis equal;
axis([x0(1)-1, x(1,end)+1, 0, x0(2)+1]);

% Start the animated GIF.
animation_name = fullfile('..', 'jade', 'img', 'ekf_demo_animation.gif');
[A, map] = rgb2ind(frame2im(getframe()), 256);
imwrite(A, map, animation_name, 'gif', 'LoopCount', inf, 'DelayTime', 2);

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
    [A, map] = rgb2ind(frame2im(getframe()), 256);
    delay = dt;
    if k == length(t)
        delay = 2;
    end
    imwrite(A, map, animation_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    
    % Wait until it's time for the next frame.
    while toc() < t(k)
        pause(0.01);
    end
    
end

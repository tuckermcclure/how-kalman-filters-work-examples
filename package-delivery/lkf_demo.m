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
vt = sqrt(9.81 / cd);
xn = [0; 0; 0; -vt];

% Initial conditions
x0  = [-1; 8; 1; 0];
P0  = bdiag(R, 1^2 * eye(2));
z0  = x0(1:2) + mnddraw(R);
xh0 = [z0; x0(3:4) + mnddraw(P0(3:4,3:4))];

% Create the filter options.
A  = [zeros(2), eye(2); zeros(2), -cd/m*vt * eye(2)];
% B  = [zeros(2); 0*eye(2)];
% Ab = [A B; zeros(2, 6)];
% Fb = expm(Ab*dt);
% F  = Fb(1:4,1:4);
% Fu = Fb(1:4, 5:6);
F  = expm(A*dt);
u  = f(xn, [0; 0]) - xn;
H  = [eye(2), zeros(2)];
% ff = @(~, ~, x, u) F * x + Fu * u;
ff = @(~, ~, x, u) F * x + u;
hf = @(~, x, ~) H * x;
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
    [xh(:,k), P(:,:,k)] = kff(xh(:,k-1) - xn, P(:,:,k-1), u, z(:,k), options);
    xh(:,k) = xh(:,k) + xn;
end

% Tools to draw ellipses
circle  = feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));
ellipse = @(P,x) bsxfun(@plus, 3 * chol(P(1:2,1:2), 'lower') * circle, x(1:2));

% Record the results with names for comparison to the EKF.
x_lkf = x; xh_lkf = xh;

%% Final plots

% % Draw the truth, measurements, and estimates as line plots.
% clf(figure(1));
% for k = 1:4
%     subplot(4, 1, k);
%     args = {t, x(k,:), ...
%             t, xh(k,:), ':'};
%     if k <= 2
%         args = [args, {t, z(k,:), '.'}]; %#ok<AGROW>
%     end
%     plot(args{:});
% end

% Show the trajectory plot.
ell0 = ellipse(P0, xh0);
set(clf(figure(2)), 'Color', [1 1 1]);
ht = plot(x(1,:),    x(2,:), ...
          z(1,:),    z(2,:),  '.', ...
          xh(1,:),   xh(2,:), 'o', ...
          ell0(1,:), ell0(2,:));
set(ht(3:4), 'Color', 0.75 * [1 1 1]);
axis equal;
axis([x0(1)-1, x(1,end)+1, 0, x0(2)+1]);

% % Animate the descent.
% tic();
% for k = 1:length(t)
%     ell = ellipse(P(:,:,k), xh(:,k));
%     set(ht(1), 'XData', x(1,1:k),  'YData', x(2,1:k));
%     set(ht(2), 'XData', z(1,1:k),  'YData', z(2,1:k));
%     set(ht(3), 'XData', xh(1,1:k), 'YData', xh(2,1:k));
%     set(ht(4), 'XData', ell(1,:),  'YData', ell(2,:));
%     while toc() < t(k)
%         pause(0.01);
%     end
% end

% Skip the animation; it looks the same.
ell = ellipse(P(:,:,end), xh(:,end));
set(ht(4), 'XData', ell(1,:),  'YData', ell(2,:));

%% Comparison

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
%     if k <= 2
%         args = [args, {t, z(k,:) - x(k,:), '.'}]; %#ok<AGROW>
%     end
    plot(args{:});
    ylabel(sprintf('State %d', k));
end
xlabel('Time (s)');
subplot(4, 1, 1);
title('Differences Between EKF and LKF');
% legend('EKF - LKF', 'Meas. - Truth');

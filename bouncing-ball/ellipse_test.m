P = randcov(2);
[u,s] = svd(P);
C = chol(P, 'lower');
% [v,e] = eig(P)

n = 1000;
X = mnddraw(P, n);

% ell = u * sqrt(s) * feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));
ell = C * feval(@(a) [cos(a); sin(a)], linspace(0, 2*pi, 100));

% X = u * sqrt(s) * Xn 
% u * X = sqrt(s) * Xn
% sqrt(s) \ u.' * X = Xn
% 
% P = 1/n * X * X.'
%   = 1/n * u * sqrt(s) * Xn * Xn.' * sqrt(s).' * u.'
%   = u * sqrt(s) * sqrt(s).' * u.'
%   = u * s * u.'

% Xn = sqrt(s) \ u.' * X;
Xn = C \ X;
a  = sum(Xn.^2, 1) < 1;

figure(3);
clf();
plot(ell(1,:), ell(2,:), ...
     X(1,~a),  X(2,~a), 'b.', ...
     X(1,a),   X(2,a),  'r.');
axis equal;

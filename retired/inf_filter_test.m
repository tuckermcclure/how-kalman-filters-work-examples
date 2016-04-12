%%
clc;

nx = 3;
nz = 5;
n = 10;
R = diag(1:nz);
invR = inv(R);
H = randn(nz, nx);
x = randn(nx, 1);
r  = covdraw(R, n);
z = H * x * ones(1, n) + r;

% zt = L * z;
% Rt = eye(nz) = L * R * L.'
% inv(L) * inv(L).' = R
% C * C.' = R
C = chol(R, 'lower');
L = inv(C);
zt = L * z;
Ht = L * H;
Rt = L * R * L.';

x

%%
H2 = repmat(H, n, 1);
z2 = H2 * x + r(:);
invR2 = diag(repmat(diag(invR), n, 1));
x2 = (H2.' * invR2 * H2) \ (H2.' * invR2) * z2

%%
(n * (H.' * invR * H)) \ (H.' * invR) * (sum(z, 2))

%%
Y = zeros(nx);
r = zeros(nx, 1);
for k = 1:n
    Y = Y + H.' * invR * H;
    r = r + H.' * invR * z(:,k);
end
Y \ r

%%
(n * (Ht.' * Ht)) \ Ht.' * (sum(zt, 2))

%%
(n * (H.' * L.' * L * H)) \ (H.' * L.') * (sum(zt, 2))

%%
Y = zeros(nx);
r = zeros(nx, 1);
for k = 1:n
    Y = Y + Ht.' * Ht;
    r = r + Ht.' * zt(:,k);
end
Y \ r

%%
Y = zeros(nx);
r = zeros(nx, 1);
for k = 1:n
    Y = Y + H.' * L.' * L * H;
    r = r + H.' * L.' * L * z(:,k);
end
Y \ r

%%
n  = 10000;
% F  = randn(2, 1);
% R  = F * rand() * F.';
R = randcov(2, 0, 100);
dz = mnddraw(R, n);

nm = 100;
[dzx, dzy] = meshgrid(linspace(-3, 3, nm), linspace(-3, 3, nm));
pdf = mndpdf([dzx(:).'; dzy(:).'], R);
pdf = reshape(pdf, nm, nm);

clf(figure(1));
histogram2(dz(1,:), dz(2,:), 'Normalization', 'pdf');
hold on;
surf(dzx, dzy, pdf);
hold off;

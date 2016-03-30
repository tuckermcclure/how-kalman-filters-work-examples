% %%
% rng(7);
% 
% x = randn(2, 1);
% P = randcov(2);
% C = chol(P, 'lower');
% a = 0.6;
% 
% S = bsxfun(@plus, a * [[0; 0], C, -C], x);
% 
% n = 256;
% X = bsxfun(@plus, mnddraw(P, n), x);
% 
% m = 100;
% [xm, ym] = meshgrid(linspace(min(X(1,:)), max(X(1,:)), m), ...
%                     linspace(min(X(2,:)), max(X(2,:)), m));
% pm = zeros(m);
% for r = 1:m
%     for c = 1:m
%         dx = [xm(r,c); ym(r,c)] - x;
%         pm(r,c) = 1/(2*pi * sqrt(det(P))) * exp(-0.5 * dx.' * inv(P) * dx); %#ok<MINV>
%     end
% end
% 
% clf(figure(1));
% contourf(xm, ym, pm, 8);
% hold on;
% h = plot(S(1,:), S(2,:), 'ko');
% hold off;
% xlabel('State 1');
% ylabel('State 2');
% colorbar();
% title('Sigma Points and Probability Density');
% xlim = get(gca(), 'XLim');
% ylim = get(gca(), 'YLim');
% % title(colorbar(), 'Probability Density');
% legend(h, 'Sigma Points');
% 
% %%
% clf(figure(1));
% plot(X(1,:), X(2,:), 'o', ...
%      S(1,:), S(2,:), 'o');
% xlabel('State 1');
% ylabel('State 2');
% title('Sigma Points vs. Particles');
% axis([xlim ylim]);
% legend('Particles', 'Sigma Points');

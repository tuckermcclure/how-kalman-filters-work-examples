cc;
% rng(1);

nw = 5;
nd = 50;

w = rand(1, nw);
w = w./sum(w);

cdf = cumsum(w);
cdf = cdf./cdf(end);

% w
w(end)

d0 = rand(1, nd);

%%

tic();
d = d0;
d = cumsum(d);
d = d./(d(end) + rand());

x = [[0 cdf]; [cdf 1]; nan(1, nw+1)];
% x = x(:);
y = [1:nw+1; 1:nw+1; nan(1, nw+1)];
plot(x(:), y(:), ...
     d,    0*d, '.');
xlabel('Random Number');
ylabel('Index');
axis([0 1 -1 nw+1]);
legend('Map', 'Draws', 'Location', 'northwest');

indices = zeros(1, nd);
wi      = 1;
for k = 1:nd
    while d(k) > cdf(wi)
        wi = wi + 1;
    end
    indices(k) = wi;
end

h = zeros(1, nw);
for k = 1:nw
    h(k) = sum(indices == k) / nd;
end
toc();

% h
h(end)

%%

tic();
d = d0;
d = sort(d);

indices = zeros(1, nd);
wi      = 1;
for k = 1:nd
    while d(k) > cdf(wi)
        wi = wi + 1;
    end
    indices(k) = wi;
end

h = zeros(1, nw);
for k = 1:nw
    h(k) = sum(indices == k) / nd;
end
toc();

% h
h(end)
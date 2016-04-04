function h = rowsubplot(p, x, y, varargin)

    np = size(y, 1);

    h = zeros(1, np);
    for k = 1:np
        subplot(np, 1, k);
        h(k) = p(x(min(k, size(x, 1)), :), y(k,:), varargin{:});
    end

end

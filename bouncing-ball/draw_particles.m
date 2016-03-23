function hso = draw_particles(X, w, xt)

    persistent hs hv ht;
    
    hso = [];
    
    % Reset.
    if nargin == 1
        hs = []; hv = []; ht = [];
        return;
    end

    % Constants
    n = size(X, 2);
    s = 0.25; % Scaling parameter for velocity arrows
    
    % Make a list of the colors.
    c = 1 - max(0.9*w/max(w).^0.75, 0.1);
    
    % Draw the balls.
    if isempty(hs)
        hs = scatter(X(1,:), X(2,:), [], c, 'o', 'filled');
    else
        set(hs, 'XData', X(1,:), 'YData', X(2,:), 'CData', c);
    end
    
    % Draw the position arrows.
    if isempty(hv)
        hv = zeros(1, n);
        for k = 1:n
            hv(k) = plot([X(1,k), X(1,k) + s * X(3,k)], ...
                         [X(2,k), X(2,k) + s * X(4,k)], ...
                         'Color', c(k) * [1 1 1]);
        end
    else
        for k = 1:n
            set(hv(k), 'XData', [X(1,k), X(1,k) + s * X(3,k)], ...
                       'YData', [X(2,k), X(2,k) + s * X(4,k)], ...
                       'Color', c(k) * [1 1 1]);
        end
    end
    
    [~, ind] = sort(w);

    % Draw/add to the trajectories.
    if nargin >= 3
        if ~isempty(ht)
            delete(ht);
            ht = [];
        end
        if ~isempty(xt)
            ht = zeros(1, n);
            for k = ind
                ht(k) = plot(xt{k}(1,:), ...
                             xt{k}(2,:), ...
                             'Color', c(k) * [1 1 1]);
            end
        end
    elseif ~isempty(ht)
        for k = 1:n
            set(ht(k), 'Color', c(k) * [1 1 1]);
        end
    end
    
    if ~isempty(ht)
        hso = ht(1);
    end

end % draw_particles

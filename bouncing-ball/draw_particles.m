function [hso, hX] = draw_particles(X, w, xt)

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
    c = 1 - max(0.95*w/max(w).^0.85, 0.05);
    
    % Draw the balls.
    if isempty(X)
        if ~isempty(hs)
            delete(hs);
            hs = [];
        end
    else
        if isempty(hs)
            hs = scatter(X(1,:), X(2,:), [], c, 'o', 'filled');
        else
            set(hs, 'XData', X(1,:), 'YData', X(2,:), 'CData', c);
        end
    end
    
%     % Draw the velocity arrows.
%     if isempty(hv)
%         hv = zeros(1, n);
%         for k = 1:n
%             hv(k) = plot([X(1,k), X(1,k) + s * X(3,k)], ...
%                          [X(2,k), X(2,k) + s * X(4,k)], ...
%                          'Color', c(k) * [1 1 1]);
%         end
%     else
%         for k = 1:n
%             set(hv(k), 'XData', [X(1,k), X(1,k) + s * X(3,k)], ...
%                        'YData', [X(2,k), X(2,k) + s * X(4,k)], ...
%                        'Color', c(k) * [1 1 1]);
%         end
%     end
    
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
    
    hX = hs;

end % draw_particles

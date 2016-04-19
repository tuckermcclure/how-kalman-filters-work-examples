function kalman_gain_demo(make_gif)

% kalman_gain_demo
%
% This function creates an animation of the correction of the covariance
% with the Kalman gain. It shows:
% 
% * the predicted state and covariance, shown as particles,
% * the predicted measurement,
% * the covariance of the predicted measurement (also particles),
% * the additional covariance of the measurement noise, yielding the
%   innovation covariance,
% * the new measurement,
% * the likelihood of each particle in observation space (just for show),
% * the likelihood of each particle in state space (also for show),
% * the Kalman correction of the predicted state,
% * and the corresponding correction of the covariance, showing it narrow
%   in on the corrected state (and, just for fun, the corresponding motion
%   of the particles in observation space).
%
% When the input argument is true, this function will create an animated
% GIF in the 'animations' folder.

% Copyright 2016 An Uncommon Lab

    % By default, don't make the GIF.
    if nargin < 1, make_gif = false; end;

    % Set a random number seed so that this animation is repeatable.
    % (Comment this line out to see different results.)
    rng(3);

    % Define the observation parts of the system.
    R = 0.25^2 * randcov(2);
    H = randn(2);

    % Define the initial conditions.
    P     = 0.5^2 * randcov(2);
    x0    = [0; 0];
    zhat0 = H * x0;

    % Calculate the Kalman gain.
    Pxy = P * H.';
    Pyy = H * P * H.' + R;
    K = Pxy / Pyy;

    % Create the true observation.
    z = H * x0 + [0.5; 0.5];%covdraw(R);

    % Create a bunch of particles to show the uncertainty.
    n  = 1000;
    X  = covdraw(P, n);
    w0 = repmat(1/n, 1, n);

    % Calculate the predicted observations.
    Z  = H * X;
    Zn = Z + covdraw(R, n);
    Dz = bsxfun(@minus, z, Zn);

    % Determine the weights based on the differences between the predicted
    % measurements (with noise) and the actual measurement.
    w = mndpdf(Dz, Pyy);
    w = w ./ sum(w); % (and normalize)

    % Calculate the corrected state.
    xc = x0 + K * (z - zhat0);

    % Calculate the corrected particles.
    Xc = X + K * Dz;
    Zc = Zn + H * K * Dz;

    % Calculate the covariance of the corrected state using two methods:
    % Kalman's and particles
    fprintf('Corrected covariance with Kalman''s equation:\n');
    Pc = P - K * H * P
    fprintf('Covariance of corrected particles:\n');
    Xct = bsxfun(@minus, Xc, xc);
    1/n * (Xct * Xct.')

    % Create a figure with the particles on the left and the predicted
    % measurements (before noise) on the right.
    clf(set(figure(1), 'Color', [1 1 1]));
    colormap(flipud(colormap('gray')));

    blue = [0 0.447 0.741];
    red  = [0.85 0.325 0.098];
    
    % Initial estimate and particles
    subplot(1, 2, 1);
    hs1 = scatter(X(1,:), X(2,:), [], w0, '.');
    hold on;
    hx0 = scatter(x0(1), x0(2), 100, 'o', 'CData', blue);
    hold off;
    clim = [-0.05 1] * max(w);
    caxis(clim);
    axis equal;
    axis(2*[-1 1 -1 1]);
    title('State Space');
    xlabel('x_1');
    ylabel('x_2');
    % legend([hx0 hs1], 'Predicted State', 'Prediction Covariance');
    
    % Observation space
    subplot(1, 2, 2);
    caxis(clim);
    axis equal;
    axis(2*[-1 1 -1 1]);
    title('Observation Space');
    xlabel('z_1');
    ylabel('z_2');
    % legend([hz0 hs2], 'Predicted Measurement', 'Prediction Covariance');

    % Start the animated GIF.
    if make_gif
        anim_name = fullfile('animations', 'kalman_gain_animation.gif');
        [A, map] = rgb2ind(frame2im(getframe(gcf())), 256);
        imwrite(A, map, anim_name, 'gif', ...
                'LoopCount', inf, ...
                'DelayTime', 2);
    end

    % A function to grab the frame data and save it in the animated gif.
    function store_frame(dt)
        if make_gif
            [A, map] = rgb2ind(frame2im(getframe(gcf())), 256);
            imwrite(A, map, anim_name, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime',  dt);
        end
    end % store_frame

    % Look at the initial plot for a moment.
    pause(1);
    
    % Initial predicted measurement
    subplot(1, 2, 2);
    hold on;
    hz0 = scatter(zhat0(1), zhat0(2), 100, 'o', 'CData', red);
    
    pause(1);
    store_frame(1);
    
    % Initial particles
    hs2 = scatter(Z(1,:), Z(2,:), [], w0, '.');
    uistack(hz0, 'top');
    hold off;
    
    pause(1);
    store_frame(1);
    
    % Animation parameters
    ns = 25;   % Number of steps to take per animation
    dt = 0.04; % Time to wait at each step [s]

    % Add the measurement noise to the predicted measurement particles to
    % get an indication of the innovation covariance.
    for k = 1:ns
        set(hs2, 'XData', Z(1,:) + k/ns*(Zn(1,:) - Z(1,:)), ...
                 'YData', Z(2,:) + k/ns*(Zn(2,:) - Z(2,:)));
        pause(dt);
        store_frame(dt);
    end
    pause(1);
    store_frame(1);

    % Add on the new measurement.
    subplot(1, 2, 2);
    hold on;
    hz = scatter(z(1), z(2), 100, 'o', 'filled', 'CData', red);
    hold off;
    pause(1);
    store_frame(1);

    % Update the weights on the predicted measurement particles.
    for k = 1:ns
        set(hs2, 'CData', w0 + k/ns*(w - w0));
        pause(dt);
        store_frame(dt);
    end
    pause(1);
    store_frame(1);

    % Update the weights on the corresponding states.
    for k = 1:ns
        set(hs1, 'CData', w0 + k/ns*(w - w0));
        pause(dt);
        store_frame(dt);
    end
    pause(1);
    store_frame(1);

    % Show the corrected state.
    subplot(1, 2, 1);
    hold on;
    hxc = scatter(xc(1), xc(2), 100, 'o', 'filled', 'CData', blue);
    hold off;
    pause(1);
    store_frame(1);

    % Move the particles showing the covariance of the state and innovation
    % to the corrected positions.
    for k = 1:ns
        set(hs1, 'XData', X(1,:) + k/ns*(Xc(1,:) - X(1,:)), ...
                 'YData', X(2,:) + k/ns*(Xc(2,:) - X(2,:)));
        set(hs2, 'XData', Zn(1,:) + k/ns*(Zc(1,:) - Zn(1,:)), ...
                 'YData', Zn(2,:) + k/ns*(Zc(2,:) - Zn(2,:)));
        pause(dt);
        store_frame(dt);
    end
    pause(1);
    store_frame(1);

    % Add on the legend at the end.
    subplot(1, 2, 1);
    legend([hx0 hs1 hxc], ...
           'Predicted State', ...
           'Corrected Covariance', ...
           'Corrected State');
    subplot(1, 2, 2);
    legend([hz0 hs2 hz], ...
           'Predicted Measurement', ...
           'Corrected Covariance', ...
           'New Measurement');
    store_frame(2);

end % kalman_gain_demo

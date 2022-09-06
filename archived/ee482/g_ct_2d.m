function g_ct_2d (funz, t, rangz)
%{
    syms t;
    f_1 = sin(2 * t);
    f_2 = sin(3 * t);
    start = 0;
    stop = 3*pi;
    g_ct_2d([f_1, f_2], t, [start, stop; start-pi, stop+pi]);
%}
    posz = [20, 20, 700, 700];    
    bufX = .1;
    negX = -100;
    posX = 100;
    bufY = .1; 
    negY = -100;
    posY = 100;
    blips = 500;
    
    burner = size(funz);
    funs = burner(2);
    x_vals = {};
    y_vals = {};
    for ii = 1:1:funs
        x_vals{1, ii} = linspace(rangz(ii, 1), rangz(ii, 2), blips);
        y_vals{1, ii} = double(subs(funz(ii), t, x_vals{1, ii}));
    end
    
    max_x = -99999999;
    max_y = max_x;
    min_x = -1 * max_x;
    min_y = min_x;
    for ii = 1:1:funs
        x_temp = x_vals{1, ii};
        y_temp = y_vals{1, ii};
        x_temp_max = max(x_vals{1, ii}, [], 'all');
        if x_temp_max > max_x
            max_x = x_temp_max;
            if max_x > posX
                max_x = posX;
            end
        end
        x_temp_min = min(x_temp, [], 'all');
        if x_temp_min < min_x
            min_x = x_temp_min;
            if min_x < negX
                min_x = negX;
            end
        end
        y_temp_max = max(y_temp, [], 'all');
        if y_temp_max > max_y
            max_y = y_temp_max;
            if max_y > posY
                max_y = posY;
            end
        end
        y_temp_min = min(y_temp, [], 'all');
        if y_temp_min < min_y
            min_y = y_temp_min;
            if min_y < negY
                min_y = negY;
            end
        end
    end
    
    x_rng = [min_x - bufX, max_x + bufX];
    y_rng = [min_y - bufY, max_y + bufY];
    x_ax = linspace(x_rng(1), x_rng(2), blips);
    y_ax = linspace(y_rng(1), y_rng(2), blips);
    
    figure('Position', posz);
    hold on;
    grid on;
    view(2);
    tiStr = "your function(s)";
    title(tiStr, 'fontsize', 26); 
    xlabel('x', 'fontsize', 18);
    ylabel('Y = f(x)', 'fontsize', 18);   
    xlim([x_rng(1), x_rng(2)]);
    ylim([y_rng(1), y_rng(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);        
    plot(x_rng(2), 0      , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0      , y_rng(2), 'y.', 'markersize', 20, 'linewidth', 10); 
    
    for ii = 1:1:funs
        if ii == 1
            plot(x_vals{1, ii}, y_vals{1, ii}, 'r','LineWidth', 3);
        elseif ii == 2
            plot(x_vals{1, ii}, y_vals{1, ii}, 'c','LineWidth', 3);
        elseif ii == 3
            plot(x_vals{1, ii}, y_vals{1, ii}, 'g','LineWidth', 3);
        else
            plot(x_vals{1, ii}, y_vals{1, ii}, 'b','LineWidth', 3);
        end
    end
    
    % plot custom points here if needed
end

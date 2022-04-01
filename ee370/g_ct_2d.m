function g_ct_2d (x_val, y_vals)
%{
    used to graph a continous 2D function
    this function does not return a value
    can handle multipe functions, keep x_val ==
    control big spikes if they show up
        
        x_val = [start, ..., stop]  "x"                1 x N 
        y_rng = [val_0, ..., val_n ; ... ; ...] "y"    M x N
        !!! equal columns for all vectors !!!

    ie:
        syms s;
        start = -6;
        stop = 6;
        dots = 100;
        push_s = linspace(start, stop, dots);

        Fs = cos(s);
        push_F = double(subs(Fs, s, push_s));
        Gs = sin(s);
        push_G = double(subs(Gs, s, push_s));
        Hs = Fs * Gs * s^2;
        push_H = double(subs(Hs, s, push_s));
        g_ct_2d(push_s, [push_F; push_G; push_H]); 
%}

    posz = [20, 20, 700, 700];    
    bufX = .1;
    negX = -10;
    posX = 10;
    bufY = .1; 
    negY = -10;
    posY = 10;
    blips = 100;
    x_rng = [negX - bufX, posX + bufX];
    y_rng = [negY - bufY, posY + bufY];
    
    if min(x_val, [], 'all') > x_rng(1)
        x_rng(1) = min(x_val, [], 'all') - bufX;
    end
    if max(x_val, [], 'all') < x_rng(2)
        x_rng(2) = max(x_val, [], 'all') + bufX;
    end
    if min(y_vals, [], 'all') > y_rng(1)
        y_rng(1) = min(y_vals, [], 'all') - bufY;
    end
    if max(y_vals, [], 'all') < y_rng(2)
        y_rng(2) = max(y_vals, [], 'all') + bufY;
    end
    
    x_ax = linspace(x_rng(1), x_rng(2), blips);
    y_ax = linspace(y_rng(1), y_rng(2), blips);
    [rows, cols] = size(y_vals)
    
    if size(x_val) ~= size(y_vals)
        fprintf("\n\tmismatched function\n");
        return;
    end
    if rows < 1
        fprintf("\nno data to plot\n");
        return;
    end
    if cols < 2
        fprintf("\nno data to plot\n");
        return;
    end
    
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
    
    for ii = 1:1:rows
        if ii == 1
            plot(x_val, y_vals(ii,:), 'r','LineWidth', 3);
        elseif ii == 2
            plot(x_val, y_vals(ii,:), 'c','LineWidth', 3);
        elseif ii == 3
            plot(x_val, y_vals(ii,:), 'g','LineWidth', 3);
        else
            plot(x_val, y_vals(ii,:), 'b','LineWidth', 3);
        end
    end
    
    % plot custom points here if needed
    
end

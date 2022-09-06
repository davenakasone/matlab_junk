function fun_image_z(tran, zformat, rng_rx, rng_ty)
%{
    "tran" is a transformation rule  ie  f(z) = z^2     ensure that f'(z_0) not 0 for range!!!
    "zformat" z=r e^jt  or z=x+jy  ...determines operations   1:rec 2:polar
    "rng_rx"   could be "r" or "x" range
    "rng_ty"   could be "theta" or "y" range
%}
    global Z;
    
    buf = 1;
    pos_z = [20, 20, 600, 600];
    pos_w= [620, 20, 600, 600];
    shade = 512;
    fill = 512;
    rng_x = [999,-999];
    rng_y = rng_x;
    rng_u = rng_x;
    rng_v = rng_x;
    
    if zformat == 1
        if rng_rx(1) == rng_rx(2)
            val_rx = zeros(1,shade) + rng_rx(1);
            tit_rx = sprintf('%x = .2f',rng_rx(1));
        elseif rng_rx(1) < rng_rx(2)
            val_rx = linspace(rng_rx(1), rng_rx(2), shade);
            tit_rx = sprintf('%.2f  <  x  <  %.2f', rng_rx(1), rng_rx(2));
        else
            val_rx = linspace(rng_rx(2), rng_rx(1), shade);
            tit_rx = sprintf('%.2f  <  x  <  %.2f', rng_rx(2), rng_rx(1));
        end
        if rng_ty(1) == rng_ty(2)
            val_ty = zeros(shade) + rng_ty(1);
            tit_ty = sprintf('y = %s',rng_ty(1));
        elseif rng_ty(1) < rng_ty(2)
            val_ty = linspace(rng_ty(1), rng_ty(2), shade);
            tit_ty = sprintf('%.2f  <  y  <  %.2f', rng_ty(1), rng_ty(2));
        else
            val_ty = linspace(rng_ty(2), rng_ty(1), shade);
            tit_ty = sprintf('%.2f  <  y  <  %.2f', rng_ty(2), rng_ty(1));
        end
        tit_z = sprintf(' z = x + jy  ::  %s  ::  %s', tit_rx, tit_ty);
        tit_w = sprintf('w = f(z) = %s', tran);
    end
    
    if zformat == 2
        if rng_rx(1) == rng_rx(2)
            val_rx = zeros(1,shade) + rng_rx(1);
            tit_rx = sprintf('r = %.2f',rng_rx(1));
        elseif rng_rx(1) < rng_rx(2)
            val_rx = linspace(rng_rx(1), rng_rx(2), shade);
            tit_rx = sprintf('%.2f  <  r  <  %.2f', rng_rx(1), rng_rx(2));
        else
            val_rx = linspace(rng_rx(2), rng_rx(1), shade);
            tit_rx = sprintf('%.2f  <  r  <  %.2f', rng_rx(2), rng_rx(1));
        end
        if rng_ty(1) == rng_ty(2)
            val_ty = zeros(shade) + rng_ty(1);
            tit_ty = sprintf('th = %s',rng_ty(1));
        elseif rng_ty(1) < rng_ty(2)
            val_ty = linspace(rng_ty(1), rng_ty(2), shade);
            tit_ty = sprintf('%.2f  <  th  <  %.2f', rng_ty(1), rng_ty(2));
        else
            val_ty = linspace(rng_ty(2), rng_ty(1), shade);
            tit_ty = sprintf('%.2f  <  th  <  %.2f', rng_ty(2), rng_ty(1));
        end
        tit_z = sprintf(' z = r exp(jt)  ::  %s  ::  %s ', tit_rx, tit_ty);
        tit_w = sprintf('w = f(z) = %s', tran);
    end
    
    [val_rx_rx, val_ty_ty] = meshgrid(val_rx, val_ty);
    map_x = zeros(shade, shade);
    map_y = map_x;
    map_u = map_x;
    map_v = map_x;
    for row = 1:shade
        for col = 1:shade
            if zformat == 1
                temp_z = val_rx_rx(row,col) + 1j * val_ty_ty(row,col);
            end
            if zformat == 2
                temp_z = val_rx_rx(row,col) * exp( 1j * val_ty_ty(row,col) );
            end
            temp_w = subs(tran, Z, temp_z);
            map_x(row,col) = double(real(temp_z));
            map_y(row,col) = double(imag(temp_z));
            map_u(row,col) = double(real(temp_w));
            map_v(row,col) = double(imag(temp_w));
            if map_x(row,col) < rng_x(1)
                rng_x(1)= map_x(row,col);
            end
            if map_x(row,col) > rng_x(2)
                rng_x(2) = map_x(row,col);
            end
            if map_y(row,col) < rng_y(1)
                rng_y(1)= map_y(row,col);
            end
            if map_y(row,col) > rng_y(2)
                rng_y(2) = map_y(row,col);
            end
            if map_u(row,col) < rng_u(1)
                rng_u(1)= map_u(row,col);
            end
            if map_u(row,col) > rng_u(2)
                rng_u(2) = map_u(row,col);
            end
            if map_v(row,col) < rng_v(1)
                rng_v(1)= map_v(row,col);
            end
            if map_v(row,col) > rng_v(2)
                rng_v(2) = map_v(row,col);
            end
        end
    end
    %rng_x = [ min(map_x, [], 'all'), max(map_x, [], 'all') ];
    %rng_y = [ min(map_y, [], 'all'), max(map_y, [], 'all') ];
    %rng_u = [ min(map_u, [], 'all'), max(map_u, [], 'all') ];
    %rng_v = [ min(map_v, [], 'all'), max(map_v, [], 'all') ];
    ax_x = linspace(rng_x(1)-buf, rng_x(2)+buf, fill);
    ax_y = linspace(rng_y(1)-buf, rng_y(2)+buf, fill);
    ax_u = linspace(rng_u(1)-buf, rng_u(2)+buf, fill);
    ax_v = linspace(rng_v(1)-buf, rng_v(2)+buf, fill);
    
    figure('Name', 'z inputs',...
               'Position', pos_z,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    view(2);
    title(tit_z, 'FontSize', 18);
    xlabel('x = real(z)', 'FontSize', 14);
    ylabel('y = imag(z)', 'FontSize', 14);
    xlim([double(rng_x(1)-buf), double(rng_x(2))+buf]);
    ylim([double(rng_y(1))-buf, double(rng_y(2))+buf]);
    plot(ax_x   , 0*ax_x , 'k' , 'linewidth' , 1);
    plot(0*ax_y , ax_y   , 'k' , 'linewidth' , 1);        
    plot(rng_x(2)+buf, 0             , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0            , rng_y(2)+buf , 'y.', 'markersize', 20, 'linewidth', 10); 
    text(rng_x(2)+buf, 0           , 'x', 'FontSize', 16);
    text(0           , rng_y(2)+buf, 'y', 'FontSize', 16);
    plot(map_x, map_y, 'r.', 'MarkerSize', 10, 'LineWidth', 5);
    hold off;
    
    figure('Name', 'w outputs',...
               'Position', pos_w,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    view(2);
    title(tit_w, 'FontSize', 18);
    xlabel('u = real(w)', 'FontSize', 14);
    ylabel('v = imag(w)', 'FontSize', 14);
    xlim([double(rng_u(1)-buf), double(rng_u(2)+buf)]);
    ylim([double(rng_v(1))-buf, double(rng_v(2))+buf]);
    plot(ax_u   , 0*ax_u , 'k' , 'linewidth' , 1);
    plot(0*ax_v , ax_v   , 'k' , 'linewidth' , 1);        
    plot(rng_u(2)+buf, 0             , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0            , rng_v(2)+buf , 'y.', 'markersize', 20, 'linewidth', 10); 
    text(rng_u(2)+buf, 0           , 'u', 'FontSize', 16);
    text(0           , rng_v(2)+buf, 'v', 'FontSize', 16);
    plot(map_u, map_v, 'b.', 'MarkerSize', 10, 'LineWidth', 5);
end


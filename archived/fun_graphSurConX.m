function fun_graphSurConX(funs, select)
%{
    [funs] = [ funz1, funz2, ..., funzX]     how ever many you want
    select = 1, just graph real part                         RED
    select = 2, just graph imaginary part                    BLUE
    select = 3, just graph abs/mag                           GREEN
    select = 4, just graph phase/angle ...the principal      '#EDB120' orange...
    select = 5, graph real, imag, abs
    select = 6, graph real, imag, abs, and phase 

    be sure to subs all functions coming in to take global Zxy  so Z = X + jY
%}
buf = 1;         % CHANGE
global X; 
global Y;
bump = 30;
posS = [30, 0, 700, 700];
posC = [posS(1)+posS(3)+bump*3, 100, 600, 600];
temp = size(funs);
siz = temp(2);
xRng = [-7, 7];   % CHANGE CHANGE CHANGE    range x inputs take [min, max]
yRng = [-3, 3];   % CHANGE CHANGE CHANGE    range y inputs take [min, max]
points = 25;      % CHANGE CHANGE CHANGE    meshgrid( points, points, points)

for idx = siz:-1:1
    posSt = posS;
    posCt = posC;
    posSt(2) = posSt(2)+idx*bump;
    posCt(2) = posCt(2)+idx*bump;
    tis = sprintf('f(z) = %s', funs(1,idx));
    
    fluffx = linspace(xRng(1)-buf, xRng(2)+buf, points);
    fluffy = linspace(yRng(1)-buf, yRng(2)+buf, points);
    xCord = linspace(xRng(1), xRng(2), points);
    yCord = linspace(yRng(1), yRng(2), points);
    [x, y] = meshgrid(xCord, yCord);
    funz = zeros(points, points);
    riderZ = zeros(1,points);
    
    for m = 1:points
        riderZ(1,m) = subs(funs(1,idx), [X, Y], [xCord(1,m), 0]);
        for n = 1:points
            funz(m,n) = subs(funs(1,idx), [X, Y], [x(m,n), y(m,n)]);  % populate target function
            
        end
    end
 
    funz_r = real(funz);  
    funz_i = imag(funz);
    funz_m = abs(funz);
    funz_p = angle(funz);
    
    figure('Name', sprintf('function %d', idx),...
           'Position', posSt,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    view(125,30); % CHANGE
    title(tis, 'FontSize', 16);
    xlabel('real  x  input');
    ylabel('imag  y  input'); 
    xlim([xRng(1)-buf, xRng(2)+buf]);
    ylim([yRng(1)-buf, yRng(2)+buf]);
    plot3(fluffx  , 0*fluffx, 0*fluffx, 'k', 'linewidth', 1);
    plot3(0*fluffy, fluffy  , 0*fluffy, 'k', 'linewidth', 1);      
    plot3(xRng(2)+buf , 0          , 0 ,'y.', 'markersize', 20, 'linewidth', 10); % +x   real 
    plot3(0           , yRng(2)+buf, 0 ,'y.', 'markersize', 20, 'linewidth', 10); % +y   imag 
    text(xRng(2)+buf, 0          , 0            , '+ X');
    text(0          , yRng(2)+buf, 0            , '+ Y');
    
    if select == 1
        zRng = [ min( [min(funz_r, [], 'all'), min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [max(funz_r, [], 'all'), max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('real output');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / real 
        
        plot3(xCord, 0*xCord, riderZ, 'b--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_r,...
            'FaceColor', 'r',...
            'FaceAlpha', .6,...
            'EdgeColor', 'none'); % red = real
        mesh(x, y, funz_r,...
            'FaceColor', 'none',...
            'EdgeColor', 'w', 'FaceAlpha', .9,...
            'LineStyle',':',...
            'LineWidth',.5);
    end
    if select == 2
        zRng = [ min( [min(funz_i, [], 'all'), min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [max(funz_i, [], 'all'), max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('imag output');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / imag 
        
        plot3(xCord, 0*xCord, riderZ, 'r--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_i,...
            'FaceColor', 'b',...
            'FaceAlpha', .6,...
            'EdgeColor', 'none'); % blue = imag
        mesh(x, y, funz_i,...
            'FaceColor', 'none',...
            'EdgeColor', 'w',...
            'FaceAlpha', .9,...
            'LineStyle',':',...
            'LineWidth',.5);
    end
    if select == 3
        zRng = [ min( [min(funz_m, [], 'all'), min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [max(funz_m, [], 'all'), max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('abs output');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / abs 
        
        plot3(xCord, 0*xCord, riderZ, 'r--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_m,...
            'FaceColor', 'g',...
            'FaceAlpha', .6,...
            'EdgeColor', 'none'); % green = abs/mag
        mesh(x, y, funz_m,...
            'FaceColor', 'none',...
            'EdgeColor', 'w',...
            'FaceAlpha', 1,...
            'LineStyle',':',...
            'LineWidth',1);
    end
    if select == 4
        zRng = [ min( [min(funz_p, [], 'all'), min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [max(funz_p, [], 'all'), max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('phase angle');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / in radians
        
        plot3(xCord, 0*xCord, riderZ, 'b--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_p,...
            'FaceColor', '#EDB120',...
            'FaceAlpha', .6,...
            'EdgeColor', 'none'); % orange = phase
        mesh(x, y, funz_p,...
            'FaceColor', 'none',...
            'EdgeColor', 'w',...
            'FaceAlpha', .9,...
            'LineStyle',':',...
            'LineWidth',1);
    end
    if select == 5
        zRng = [ min( [ min(funz_r, [], 'all'),...
                       min(funz_i, [], 'all'),...
                       min(funz_m, [], 'all'),...
                       min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [ max(funz_r, [], 'all'),...
                       max(funz_i, [], 'all'),...
                       max(funz_m, [], 'all'),...
                       max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('real, imag, abs');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / real, imag, mag
        
        plot3(xCord, 0*xCord, riderZ, 'c--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_r, 'FaceColor', 'r', 'FaceAlpha', .6, 'EdgeColor', 'none'); % red = real    
        surf(x, y, funz_i, 'FaceColor', 'b', 'FaceAlpha', .6, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, funz_m, 'FaceColor', 'g', 'FaceAlpha', .6, 'EdgeColor', 'none'); % green = abs/mag
        mesh(x, y, funz_r, 'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, funz_i, 'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, funz_m, 'FaceColor', 'g','EdgeColor', 'g', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
    end
    if select == 6
        zRng = [ min( [ min(funz_r, [], 'all'),...
                       min(funz_i, [], 'all'),...
                       min(funz_m, [], 'all'),...
                       min(funz_p, [], 'all'),...
                       min(riderZ, [], 'all') ], [], 'all')-buf,...
                 max( [ max(funz_r, [], 'all'),...
                       max(funz_i, [], 'all'),...
                       max(funz_m, [], 'all'),...
                       max(funz_p, [], 'all'),...
                       max(riderZ, [], 'all') ], [], 'all')+buf];
        zax = linspace(zRng(1), zRng(2), points);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        zlabel('real, imag, abs, phase');
        zlim([zRng(1), zRng(2)]);
        plot3(0 , 0 , zRng(2),'y.', 'markersize', 20, 'linewidth', 10); % z / real, imag, mag, ang
        
        plot3(xCord, 0*xCord, riderZ, 'c--', 'LineWidth', 2); % the actual function, in real
        surf(x, y, funz_r, 'FaceColor', 'r', 'FaceAlpha', .6, 'EdgeColor', 'none'); % red = real    
        surf(x, y, funz_i, 'FaceColor', 'b', 'FaceAlpha', .6, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, funz_m, 'FaceColor', 'g', 'FaceAlpha', .6, 'EdgeColor', 'none'); % green = abs/mag
        surf(x, y, funz_p, 'FaceColor', '#EDB120', 'FaceAlpha', .6, 'EdgeColor', 'none'); % orange = phase
        mesh(x, y, funz_r, 'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, funz_i, 'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, funz_m, 'FaceColor', 'g','EdgeColor', 'g', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, funz_p, 'FaceColor', '#EDB120','EdgeColor', '#EDB120', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
    end
    hold off;
        
    
    figure('Name', sprintf('contours %d', idx),...
           'Position', posCt,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    axis equal;
    view(2); % above
    title(tis,'FontSize', 16);
    xlabel('real  x');
    ylabel('imag  y'); 
    xlim([xRng(1)-buf, xRng(2)+buf]);
    ylim([yRng(1)-buf, yRng(2)+buf]);
    plot(fluffx  , 0*fluffx, 'k', 'linewidth', 1);
    plot(0*fluffy, fluffy  , 'k', 'linewidth', 1);
    plot(xRng(2)+buf, 0          , 'y.', 'markersize', 20, 'linewidth', 10); % +x / real 
    plot(0          , yRng(2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / imag 
    % CHANGE CHANGE CHANGE
    %[pR, pI] = contour(x, y, f, [-3, -2, -1, 0, 1, 2, 3]);  select specified levels 
                            % (x, y, f, 5)   to auto select 5 levels
                            % CHANGE CHANGE CHANGE
    [fxR, fyR] = contour(x, y, funz_r, [-5, -1, 0, 1, 5],'r', 'LineWidth', 2); % real contour
    clabel(fxR, fyR);
    [fxI, fyI] = contour(x, y, funz_i, [-5, -1, 0, 1, 5],'b', 'LineWidth', 2); % imag contour
    clabel(fxI, fyI);
    hold off;
end
end

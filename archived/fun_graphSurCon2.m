function fun_graphSurCon2(funIn1, funIn2, select)
%{
    funIn1 is a surface, funIn2 is a mesh/grid
    select = 1, just graph real part                         RED
    select = 2, just graph imaginary part                    BLUE
    select = 3, just graph abs/mag                           GREEN
    select = 4, just graph phase/angle ...the principal      '#EDB120' orange...
    select = 5, graph real, imag, abs
    select = 6, graph real, imag, abs, and phase 
%}

titleStr = sprintf('%s   и   %s',funIn1, funIn2);
global X; % need X and Y as the symbols to z = x + jy
global Y;
%global Z;
    maxP = 10;       % CHANGE CHANGE CHANGE
    points = 25;   % CHANGE CHANGE CHANGE
    zf = 2;         % CHANGE CHANGE CHANGE
    xCord = linspace(-maxP, maxP, points);
    yCord = linspace(-maxP, maxP, points);
    [x, y] = meshgrid(xCord, yCord);
    f1 = zeros(points, points);
    f2 = zeros(points, points);
    for i = 1:points
        for k = 1:points
            f1(i,k) = subs(funIn1, [X, Y], [x(i,k), y(i,k)]);  % populate target function1
            f2(i,k) = subs(funIn2, [X, Y], [x(i,k), y(i,k)]);  % populate target function2
        end
    end
    
    fr1 = real(f1); % can only plot 3 vars ...2 from x, y   thrid as mag, real, imag, phase, ect
    fi1 = imag(f1);
    fmag1 = abs(f1);
    fp1 = angle(f1);
    fr2 = real(f2); 
    fi2 = imag(f2);
    fmag2 = abs(f2);
    fp2 = angle(f2);
    
    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    view(125,30);
    title(titleStr, 'fontsize', 16); % CHANGE CHANGE CHANGE
    xlabel('real');
    ylabel('imag'); 
    %zlabel('mag, real, imag....');   % CHANGE CHANGE CHANGE
    xlim([-maxP-2, maxP+2]);
    ylim([-maxP-2, maxP+2]);
    zlim([-maxP*zf-2*zf, maxP*zf+2*zf]);
    xax = linspace(-maxP-2   , maxP+2   , points);
    yax = linspace(-maxP-2   , maxP+2   , points);
    zax = linspace(-maxP*zf-2*zf, maxP*zf+2*zf, points);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(maxP-.5, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , maxP-.5, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0      , 0      , zf*maxP-.5,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real value
    % CHANGEs  ...mesh, surf, ect
    if select == 1
        surf(x, y, fr1, 'FaceColor', 'r', 'FaceAlpha', .6, 'EdgeColor', 'none'); % red = real
        mesh(x, y, fr2, 'FaceColor', 'none','EdgeColor', 'w', 'FaceAlpha', .9, 'LineStyle',':', 'LineWidth',2);
    end
    if select == 2
        surf(x, y, fi1, 'FaceColor', 'b', 'FaceAlpha', .6, 'EdgeColor', 'none'); % blue = imag
        mesh(x, y, fi2, 'FaceColor', 'none','EdgeColor', 'w', 'FaceAlpha', .9, 'LineStyle',':', 'LineWidth',2);
    end
    if select == 3
        surf(x, y, fmag1, 'FaceColor', 'g', 'FaceAlpha', .6, 'EdgeColor', 'none'); % green = abs/mag
        mesh(x, y, fmag2, 'FaceColor', [.9, .9, .9],'EdgeColor', 'w', 'FaceAlpha', .9, 'LineStyle',':', 'LineWidth',2);
    end
    if select == 4
        surf(x, y, fp1, 'FaceColor', '#EDB120', 'FaceAlpha', .6, 'EdgeColor', 'none'); % orange = phase
        mesh(x, y, fp2, 'FaceColor', [.9, .9, .9],'EdgeColor', 'w', 'FaceAlpha', .9, 'LineStyle',':', 'LineWidth',2);
    end
    if select == 5
        surf(x, y, fr1, 'FaceColor', 'r', 'FaceAlpha', .6, 'EdgeColor', 'none'); % red = real    
        surf(x, y, fi1, 'FaceColor', 'b', 'FaceAlpha', .6, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, fmag1, 'FaceColor', 'g', 'FaceAlpha', .6, 'EdgeColor', 'none'); % green = abs/mag
        mesh(x, y, fr2, 'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, fi2, 'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
        mesh(x, y, fmag2, 'FaceColor', 'g','EdgeColor', 'g', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',.1);
    end
    if select == 6
        surf(x, y, fr1, 'FaceColor', 'r', 'FaceAlpha', .6, 'EdgeColor', 'none'); % red = real    
        surf(x, y, fi1, 'FaceColor', 'b', 'FaceAlpha', .6, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, fmag1, 'FaceColor', 'g', 'FaceAlpha', .6, 'EdgeColor', 'none'); % green = abs/mag
        surf(x, y, fp1, 'FaceColor', '#EDB120', 'FaceAlpha', .6, 'EdgeColor', 'none'); % orange = phase
        mesh(x, y, fr2, 'FaceColor', 'r','EdgeColor', 'r', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',1);
        mesh(x, y, fi2, 'FaceColor', 'b','EdgeColor', 'b', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',1);
        mesh(x, y, fmag2, 'FaceColor', 'g','EdgeColor', 'g', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',1);
        mesh(x, y, fp2, 'FaceColor', '#EDB120','EdgeColor', '#EDB120', 'FaceAlpha', .2, 'LineStyle',':', 'LineWidth',1);
    end
        
    
    figure('Position',[820, 300, 500, 500]);
    hold on;
    grid on;
    axis equal;
    view(2); % above
    title('imag vs real contours ... should be othoganal', 'fontsize', 16);
    xlabel('real');
    ylabel('imag'); 
    xlim([-maxP, maxP]);
    ylim([-maxP, maxP]);
    xax = linspace(-maxP   , maxP   , points);
    yax = linspace(-maxP   , maxP   , points);
    %plot(xax  , 0*xax, 'k', 'linewidth', 1);
    %plot(0*yax, yax  , 'k', 'linewidth', 1);
    plot(maxP-.5, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot(0      , maxP-.5, 'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    % CHANGE CHANGE CHANGE
    %[pR, pI] = contour(x, y, f, [-3, -2, -1, 0, 1, 2, 3]);  select specified levels 
                            % (x, y, f, 5)   to auto select 5 levels
                            % CHANGE CHANGE CHANGE
    [frR1, frI1] = contour(x, y, fr1,5,'r'); % CHANGE CHANGE CHANGE   fr, fmag, ect
    clabel(frR1, frI1);
    [fiR2, fiI2] = contour(x, y, fi1,5,'b'); % CHANGE CHANGE CHANGE  fr, fmag, ect
    clabel(fiR2, fiI2);
    
    [frR3, frI3] = contour(x, y, fr2,5, 'r--', 'linewidth', 3); % CHANGE CHANGE CHANGE   fr, fmag, ect
    clabel(frR3, frI3);
    [fiR4, fiI4] = contour(x, y, fi2,5 ,'b--', 'linewidth', 3); % CHANGE CHANGE CHANGE  fr, fmag, ect
    clabel(fiR4, fiI4);
end

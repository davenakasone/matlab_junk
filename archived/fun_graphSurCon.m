function fun_graphSurCon(funIn, select)
%{
    select = 1, just graph real part                         RED
    select = 2, just graph imaginary part                    BLUE
    select = 3, just graph abs/mag                           GREEN
    select = 4, just graph phase/angle ...the principal      '#EDB120' orange...
    select = 5, graph real, imag, abs
    select = 6, graph real, imag, abs, and phase 
%}
    
titleStr = sprintf('%s',funIn);
global X; % need X and Y as the symbols to z = x + jy
global Y;
%global Z;
    maxP = 10;       % CHANGE CHANGE CHANGE
    points = 25;     % CHANGE CHANGE CHANGE
    zf = 2;          % CHANGE CHANGE CHANGE
    xCord = linspace(-maxP, maxP, points);
    yCord = linspace(-maxP, maxP, points);
    [x, y] = meshgrid(xCord, yCord);
    f = zeros(points, points);
    for i = 1:points
        for k = 1:points
            f(i,k) = subs(funIn, [X, Y], [x(i,k)+eps, y(i,k)+eps]);  % populate target function
        end
    end
    
    fr = real(f); % can only plot 3 vars ...2 from x, y   thrid as mag, real, imag, phase, ect
    fi = imag(f);
    fmag = abs(f);
    fp = angle(f);
    
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
    
    if select == 1
        surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none'); % red = real
    end
    if select == 2
        surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none'); % blue = imag
    end
    if select == 3
        surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .2, 'EdgeColor', 'none'); % green = abs/mag
    end
    if select == 4
        surf(x, y, fp, 'FaceColor', '#EDB120', 'FaceAlpha', .2, 'EdgeColor', 'none'); % orange phase
    end
    if select == 5
        surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none'); % red = real    
        surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .2, 'EdgeColor', 'none'); % green = abs/mag
    end
    if select == 6
        surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none'); % red = real    
        surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none'); % blue = imag
        surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .2, 'EdgeColor', 'none'); % green = abs/mag
        surf(x, y, fp, 'FaceColor', '#EDB120', 'FaceAlpha', .2, 'EdgeColor', 'none'); % orange phase
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
    plot(xax  , 0*xax, 'k', 'linewidth', 1);
    plot(0*yax, yax  , 'k', 'linewidth', 1);
    plot(maxP-.5, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot(0      , maxP-.5, 'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    % CHANGE CHANGE CHANGE
    %[pR, pI] = contour(x, y, f, [-3, -2, -1, 0, 1, 2, 3]);  select specified levels 
                            % (x, y, f, 5)   to auto select 5 levels
                            % CHANGE CHANGE CHANGE
    [frR1, frI1] = contour(x, y, fr,5,'r:', 'linewidth', 3); % CHANGE CHANGE CHANGE   fr, fmag, ect
    clabel(frR1, frI1);
    [fiR2, fiI2] = contour(x, y, fi,5,'b:', 'linewidth', 3); % CHANGE CHANGE CHANGE  fr, fmag, ect
    clabel(fiR2, fiI2);
end

   
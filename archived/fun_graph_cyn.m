function  fun_graph_cyn(axis,rad, cent, hgt, pts)
%{
   axis: 'x', 'y', 'z' ... where the major axis of cylinder is parallel 
    going to add a generalized space curve in eventually 
        that way you can wrap a cylinder around a helix or something
   
   rad:  how big   [ 1 x 1 ]

   cent:  [ 1 x 2]
            axis = 'z' ,   x = cent(1), y = cent(2)
            axis = 'y' ,   x = cent(1), z = cent(2)
            axis = 'x' ,   y = cent(1), z = cent(2)

    hgt:  [ 1 x 2 ]    [ top, bottom ]

    pts: [ n x 3 ]    points must be in rectangular   will be plotted
%}
buf = 3;
points = 20;

if axis == 'z'
    tiStr = sprintf('cyn, %s major axis, x-center: %.1f, y-center: %.1f, ht: %.1f',...
        axis, cent(1), cent(2), hgt(2)-hgt(1));
    
    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    view(125,30); % 3  also good
    title(tiStr, 'fontsize', 16); 
    xlabel('x axis');
    ylabel('y axis'); 
    zlabel('z axis');  
    xlim([cent(1)-rad-buf, cent(1)+rad+buf]);
    ylim([cent(2)-rad-buf, cent(2)+rad+buf]);
    zlim([hgt(1)-buf, hgt(2)+buf]);
    xax = linspace(cent(1)-rad-buf   , cent(1)+rad+buf   , points);
    yax = linspace(cent(2)-rad-buf   , cent(2)+rad+buf   , points);
    zax = linspace(hgt(1)-buf, hgt(2)+buf, points);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(cent(1)+rad+buf-1, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , cent(2)+rad+buf-1, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0      , 0      , hgt(2)+buf-1,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
    
    szchk = size(pts);
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'b.', 'markersize', 20);
        end
    end
    
    [x, y, z] = cylinder(rad);
    x = x + cent(1);
    y = y + cent(2);
    z1 = z * hgt(1);
    z2 = z * hgt(2);
    surf(x, y, z1, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    surf(x, y, z2, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');

    points = 10*points; % don't make 2 big
    z1 = zeros(points, points) + hgt(1);
    z2 = zeros(points, points) + hgt(2);
    %z2 = zeros(siz(1), siz(2)) + hgt(2);
    %fill3(x, y, z1,'r');
    %fill3(x(1,:), y(1,:), z1(1,:), 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    %fill3( [x(1,:), y(1,:), z2(1,:), 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    %fill3( x(2,:), y(2,:), z2(2,:), 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    t = linspace(0, 2*pi, points);
    x = rad * cos(t) + cent(1);
    y = rad * sin(t) + cent(2);
    [x, y] = meshgrid(x,y);
    for i = 1:points
        for k = 1:points
            if x(i,k)^2 + y(i,k)^2 >= rad^2
                x(i,k) = nan;
                y(i,k) = nan;
            end
        end
    end
    surf(x, y, z1, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    surf(x, y, z2, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    plot3([rad/(2^.5), rad/(2^.5)], [rad/(2^.5), rad/(2^.5)], [hgt(1), hgt(2)], 'r');
    plot3([-rad/(2^.5), -rad/(2^.5)], [rad/(2^.5), rad/(2^.5)], [hgt(1), hgt(2)], 'r');
    plot3([rad/(2^.5), rad/(2^.5)], [-rad/(2^.5), -rad/(2^.5)], [hgt(1), hgt(2)], 'r');
    plot3([-rad/(2^.5), -rad/(2^.5)], [-rad/(2^.5), -rad/(2^.5)], [hgt(1), hgt(2)], 'r');
    plot3([rad/(2^.5), -rad/(2^.5)], [rad/(2^.5), -rad/(2^.5)], [hgt(2), hgt(2)], 'r');
    plot3([-rad/(2^.5), rad/(2^.5)], [rad/(2^.5), -rad/(2^.5)], [hgt(2), hgt(2)], 'r');
    plot3([rad/(2^.5), -rad/(2^.5)], [rad/(2^.5), -rad/(2^.5)], [hgt(1), hgt(1)], 'r');
    plot3([-rad/(2^.5), rad/(2^.5)], [rad/(2^.5), -rad/(2^.5)], [hgt(1), hgt(1)], 'r');
    plot3([cent(1), cent(1)+rad/2], [cent(2),cent(2)+rad*(sqrt(3)/2)], [hgt(1), hgt(1)],...
        'g', 'linewidth', 2);
    radStr = sprintf('r = %.2f', rad);
    text(.2+cent(1)+rad/2,.2+cent(2)+rad*(sqrt(3)/2), hgt(1)-.5, radStr);
    plot3([cent(1)+rad*(sqrt(3)/2), cent(1)+rad*(sqrt(3)/2)], [cent(2)+rad/2,cent(2)+rad/2],...
        [hgt(1), hgt(2)],'g', 'linewidth', 2);
    htStr = sprintf('h = %.2f', hgt(2)-hgt(1));
    text(.2+cent(1)+rad*(sqrt(3)/2),.2+cent(2)+rad/2, hgt(2)+.5, htStr);
    
else
    fprintf('this function still has some work to be done');
end
end


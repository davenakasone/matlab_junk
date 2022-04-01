function  fun_graph_norm(funIn, pts)
%{
    funIn comes in specified by globals f(rx, ry) = rz ...improve later
        it's a 2 var function desinged to be graphed on a 3D plot

    pts (if any) will be plotted

    the function gets plotted as a blue surface,
        surfaces norms are also ploted...they are the normal vectors to the function

    you can handle different vector fields with a type var
        maybe add a center also

    or just make it a copy/paste -> graph function easy to do

%}
global rx; global ry; global rz % these are what makes the transfer happen add others as needed
buf = 1;
maxP = 5;
points = 15;
tiStr = sprintf('f(rx, ry) = %s', funIn);

    [x, y] = meshgrid( linspace(-maxP, maxP, points), linspace(-maxP, maxP, points));
    z = zeros(points, points);
    for i = 1:points
        for k = 1:points
            z(i, k) = subs(funIn, [rx, ry], [x(i, k), y(i, k)] );
        end
    end
    [u, v, w] = surfnorm(z);
    
    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    axis equal;
    view(125,30); % 3  also good
    title(tiStr, 'fontsize', 16); 
    xlabel('x axis');
    ylabel('y axis'); 
    zlabel('z axis');  
    xlim([-maxP-buf, maxP+buf]);
    ylim([-maxP-buf, maxP+buf]);
    zlim([-maxP-buf, maxP+buf]);
    xax = linspace(-maxP-buf   , maxP+buf   , points);
    yax = linspace(-maxP-buf   , maxP+buf   , points);
    zax = linspace(-maxP-buf   , maxP+buf   , points);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(maxP, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , maxP, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0      , 0   , maxP      ,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
            
    myQ = quiver3(x, y, z, u, v, w);
    myQ.ShowArrowHead = 'off';
    myQ.Marker = 'x';
    myQ.LineWidth = 2;
    
    surf(x, y, z, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none');
    szchk = size(pts);
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'g.', 'markersize', 20);
        end
    end

end







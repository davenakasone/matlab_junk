function  fun_graph_path3(pts, Vf)
%{
    pts come in specifing vertecies of graph
        row 1 => pt1 , row2 => pt 2,   ....   (x, y, z) ....add cyn and sph later

    if a vector field comes in, that will get plotted also

%}
global rx; global ry; global rz % these are what makes the transfer happen add others as needed
buf = 1;
maxPx = 0;
maxPy = 0;
maxPz = 0;
points = 5;
tiStr = sprintf('map и VF');

szComp = size(pts);
for i = 1:szComp(1)
    if maxPx < abs(pts(i,1))
        maxPx = abs(pts(i,1));
    end
    if maxPy < abs(pts(i,2))
        maxPy = abs(pts(i,2));
    end
    if maxPz < abs(pts(i,3))
        maxPz = abs(pts(i,3));
    end
end
    
    szVf = size(Vf);
    if szVf(2) > 2
    [x, y, z] = meshgrid(linspace(-maxPx, maxPx, points),...
        linspace(-maxPy, maxPy, points), linspace(-maxPz, maxPz, points));
    P = zeros(points, points, points);
    Q = zeros(points, points, points);
    R = zeros(points, points, points);
    for k = 1:points
        for m = 1:points
            for n = 1:points
                P(k, m, n) = subs(Vf(1), [rx, ry, rz], [x(k,m,n), y(k,m,n), z(k,m,n)]);
                Q(k, m, n) = subs(Vf(2), [rx, ry, rz], [x(k,m,n), y(k,m,n), z(k,m,n)]);
                R(k, m, n) = subs(Vf(3), [rx, ry, rz], [x(k,m,n), y(k,m,n), z(k,m,n)]);
            end
        end
    end
    end

    
    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    axis equal;
    view(125,30); % 3  also good
    title(tiStr, 'fontsize', 16); 
    xlabel('x axis');
    ylabel('y axis'); 
    zlabel('z axis');  
    xlim([-maxPx-buf, maxPx+buf]);
    ylim([-maxPy-buf, maxPy+buf]);
    zlim([-maxPz-buf, maxPz+buf]);
    xax = linspace(-maxPx-buf   , maxPx+buf   , points);
    yax = linspace(-maxPy-buf   , maxPy+buf   , points);
    zax = linspace(-maxPz-buf   , maxPz+buf   , points);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(maxPx+buf, 0 , 0,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0, maxPy+buf, 0,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0, 0, maxPz+buf,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
            
    myQ = quiver3(x, y, z, P, Q, R,'b');
    myQ.ShowArrowHead = 'on';
    %myQ.Marker = 'x';
    myQ.LineWidth = 1;
    
    for i = 1:szComp(1)
        if i < szComp(1)
            quiver3(pts(i,1), pts(i,2), pts(i,3),...
                pts(i+1,1)-pts(i,1), pts(i+1,2)-pts(i,2), pts(i+1,3)-pts(i,3),...
                    'r', 'linewidth',2 );
        else
            quiver3(pts(i,1), pts(i,2), pts(i,3),...
                pts(1,1)-pts(i,1), pts(1,2)-pts(i,2), pts(1,3)-pts(i,3),...
                    'r', 'linewidth',2 );
        end
    end
    
    for i = 1:szComp(1)
        plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'g.', 'markersize', 10);
    end

end


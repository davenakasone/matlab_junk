function  fun_graph_sph(rad, cent, pts)
%{
    rad [ 1 x 1 ] scalar
    cent [ 1 x 3 ]    xcord, ycord, zcord     change later if needed
    pts in rec (if any)
%}

buf = 3;
points = 20;

tiStr = sprintf('sph, radius: %.2f, x-center: %.2f, y-center: %.2f',...
        rad, cent(1), cent(2));
    
    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    axis equal;
    view(125,30); % 3  also good
    title(tiStr, 'fontsize', 16); 
    xlabel('x axis');
    ylabel('y axis'); 
    zlabel('z axis');  
    xlim([cent(1)-rad-buf, cent(1)+rad+buf]);
    ylim([cent(2)-rad-buf, cent(2)+rad+buf]);
    zlim([cent(3)-rad-buf, cent(3)+rad+buf]);
    xax = linspace(cent(1)-rad-buf   , cent(1)+rad+buf   , points);
    yax = linspace(cent(2)-rad-buf   , cent(2)+rad+buf   , points);
    zax = linspace(cent(3)-rad-buf   , cent(3)+rad+buf   , points);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(cent(1)+rad+buf-1, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , cent(2)+rad+buf-1, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0      , 0      , cent(3)+rad+buf-1,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
    
    szchk = size(pts);
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'b.', 'markersize', 20);
        end
    end
    
    [x, y, z] = sphere;
    x = x * rad + cent(1);
    y = y * rad + cent(2);
    z = z * rad + cent(3);
    surf(x, y, z, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none');
    [x, y, z] = sphere(10);
    x = x * rad + cent(1);
    y = y * rad + cent(2);
    z = z * rad + cent(3);
    [u, v, w] = surfnorm(x, y, z);
    quiver3(x, y, z, u, v, w, 0,'r');
    
    
    strl = sprintf('( %.1f, 0, 0 )', rad);
    plot3(rad, 0, 0, 'g.', 'markersize', 20);
    text(rad-.2, -.2, -.2, strl, 'FontSize', 14);
    
end


function  fun_graphEllipsoid(cent,semAx)
%{
    cent [ x coord, y cord, z cord ]
    semAx [ x axis, y axis, z axis]   the sqare root of denominator to LHS = 1
%}
figName = 'ellipsoid';

tiStr = sprintf('(x^2 - %.1f) / %.1f  +  (y^2 - %.1f) / %.1f  +  (z^2 - %.1f) / %.1f = 1',...
    cent(1), semAx(1)^2, cent(2), semAx(2)^2, cent(3), semAx(3)^2 );


buf = 1;
maxPx = abs(cent(1))+semAx(1)+buf;
maxPy = abs(cent(2))+semAx(2)+buf;
maxPz = abs(cent(3))+semAx(3)+buf;
points = 50;
    
    figure('Name', figName, 'Position',[20, 20, 700, 700]);
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
    
    plot3(cent(1), cent(2), cent(3), 'g.', 'markersize', 20);
    %{
    colormap('jet');          % ... quick way
    colorbar;
    ellipsoid(cent(1), cent(2), cent(3), semAx(1), semAx(2), semAx(3), 30); % 30 faces
    %}
    [X, Y, Z] = ellipsoid(cent(1), cent(2), cent(3), semAx(1), semAx(2), semAx(3) );
    elpS = surf(X, Y, Z, 100);
    elpS.FaceColor = 'b';
    elpS.FaceAlpha = .2;
    set(elpS, 'edgecolor', 'b');  % 'none' always good
    elpS.LineWidth = .001;
    %direction = [1, 0, 0]; % specify rotation in vector
    %rotate(elpS, direction, 45); % rotate it 45°
    
end


function  fun_graphVF(Vf, pts, levZ)
%{
    field is the rectangular vector field....going to need some work for other types
    [ rx^2 * ry, rz * 4, rx * sin(ry) ]   just an example

    pts (if any) will be plotted

    levZ is the contour height you want to view the vector feild

    you can handle different vector fields with a type var
        maybe add a center also

    or just make it a copy/paste -> graph function easy to do

%}
global rx; global ry; global rz % these are what makes the transfer happen add others as needed
buf = 2;
maxP = 5;
points = 6;
tiStr = sprintf('vector feild [ %s , %s , %s ]', Vf);
szchk = size(pts);


    % general 3D
    [x, y, z] = meshgrid(linspace(-maxP, maxP, points), linspace(-maxP, maxP, points), linspace(-maxP, maxP, points));
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
    plot3(maxP+buf, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    text(maxP+buf, 0, 0, '+ X');
    plot3(0      , maxP+buf-.2, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    text(.2, maxP+buf-.2, 0, '+ Y');
    plot3(.2      , .2   , maxP+buf      ,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
    text(.2, .2, maxP+buf-.2, '+ Z');
    quiver3(x, y, z, P, Q, R, 'r', 'AutoScale', 'on', 'AutoScaleFactor', .9);  % many options to format
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'b.', 'markersize', 20);
        end
    end
    
    
    % standard 2D      on user defined level z
    [x2d, y2d] = meshgrid(linspace(-maxP, maxP, points), linspace(-maxP, maxP, points) );
    xv = zeros(points, points);
    yv = zeros(points, points);
    for k = 1:points
        for m = 1:points
            xv(k, m) = subs(Vf(1), [rx, ry, rz], [x2d(k,m), y2d(k,m), levZ]);
            yv(k, m) = subs(Vf(2), [rx, ry, rz], [x2d(k,m), y2d(k,m), levZ]);
        end
    end
    figure();
    hold on;
    grid on;
    axis equal;
    view(2); % 2 for 2D
    tiStr = sprintf('vector feild [ %s , %s , %s ]   @  z = %.1f', Vf, levZ);
    title(tiStr, 'fontsize', 16); 
    xlabel('x axis');
    ylabel('y axis');   
    xlim([-maxP-buf, maxP+buf]);
    ylim([-maxP-buf, maxP+buf]);
    xax = linspace(-maxP-buf   , maxP+buf   , points);
    yax = linspace(-maxP-buf   , maxP+buf   , points);
    plot(xax  , 0*xax, 'k', 'linewidth', 1);
    plot(0*yax, yax  ,'k', 'linewidth', 1);        
    plot(maxP+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot(0      , maxP+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    quivH = quiver(x2d, y2d, xv, yv);
    quivH.AutoScale = 'on';
    quivH.AutoScaleFactor = .8;
    quivH.Color = 'b';
    quivH.LineWidth = 2;
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot(pts(i, 1), pts(i, 2), 'r.', 'markersize', 20);
        end
    end
    
    
    % the 3D way  for user defined level Z
    [px, py] = meshgrid(linspace(-maxP, maxP, points), linspace(-maxP, maxP, points));
    pz = zeros(points, points) + levZ;
    z = zeros(points, points, points) + levZ;
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
    figure();
    hold on;
    grid on;
    axis equal;
    view(125,30); % 3  also good
    tiStr = sprintf('vector feild [ %s , %s , %s ]   @  z = %.1f', Vf, levZ);
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
    plot3(maxP+buf, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , maxP+buf, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    plot3(0      , 0   , maxP+buf      ,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real
    quivH = quiver3(x, y, z, P, Q, R);
    quivH.AutoScale = 'on';
    quivH.AutoScaleFactor = 1;
    quivH.Color = 'b';
    quivH.LineWidth = 2;
    surf(px, py, pz, 'FaceColor', 'b', 'FaceAlpha', .1, 'EdgeColor', 'none'); % blue = imag
    if szchk(1) > 0 && szchk(2) > 2
        for i = 1:szchk(1)
            plot3(pts(i, 1), pts(i, 2), pts(i, 3), 'r.', 'markersize', 20);
        end
    end
end


%{

        it actually works fine, just easier to mesh as above
temp = linspace(maxP, -maxP, points);
x = linspace(-maxP, maxP, points);
y = linspace(-maxP, maxP, points);
[x, y] = meshgrid(x, y);
xC = zeros(points, points, points);
yC = zeros(points, points, points);
zC = zeros(points, points, points);
uC = zeros(points, points, points);
vC = zeros(points, points, points);
wC = zeros(points, points, points);
for m = 1:points
    zC(:,:,m) = temp(m);
    for n = 1:points
        for p = 1:points
            xC(m,n,p) = x(m,n);
            yC(m,n,p) = y(m,n);
        end
    end
end


figure();
hold on;
view(140,20);


                         cube of points confirmed  ...that is nice
for m = 1:points       
    for n = 1:points
        for p = 1:points
            plot3(xC(m, n, p), yC(m, n, p), zC(m, n, p), 'b.', 'markersize', 3);
        end
    end
end
%}


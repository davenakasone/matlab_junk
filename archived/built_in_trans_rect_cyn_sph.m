%{
    matlab should have some decent built in functions for spherical, cylindrical, and others
    the big 3 for 3D are spherical, cylindrical, and rectangular
    
    also going to do some 2D  ( but that is just polar)
        fyi, use retangular as your bridge

    #1 given (x, y) produce polar (theta, rho)
                                                [theta, rho] = cart2pol(x, y);
                                                [ x, y ] = pol2cart(theta, rho); 

    #2 given (x, y, z) produce cylindrical (rho, fi, z)
            note that matlab is a little different ...they use the good way (theta, rho, z)
                [theta, rho, z] = cart2pol(x, y, z);
                [x, y, z] = pol2cart(theta, rho, z);

    #3  rectangular to spherical... matlab uses (x, y, z) -> (
                azimuth lives on the xy plane...it is what fi is in the book
                but elevation is off xy plane ( -pi/2 , pi/2)    not off z  (0, pi)
                add pi/2  if less than 0
                r -> rho
                theta -> azimuth
                fi -> elevations (but adjust angle!
            [-pi , pi] on azimuth
            [pi/2 , -pi/2] on elevation
                
    
    don't forget about rad2deg and deg2rad
            
%}
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice
clc;
close all;
clearvars;
%consts = cls_CONST();
%consts.check();
%consts.help();



sel = 3;    % CHANGE to see desired calculation , see header for options



%----------------------------------------------------------------------------------------------- #1
if sel == 1
    recx = -1;         % rectangular x component                   CHANGE HERE
    recy = 1;         % rectangular y component                   CHANGE HERE
    pt = 0;           % polar theta component
    pr = 0;           % polar rho compoenent
    
    % 2D polar is easy, you can already use the regular vector functions     watch that arctan
    pt = atan(recy / recx);
    pr = sqrt( recx^2 + recy^2 );
    fprintf('the point ( %d , %d ) in rect is ( %d, %d ) in polar\n',...
        recx, recy, pt, pr);
    
    % convert polar beack to rectangular just to check
    recx = pr * cos(pt);
    recy = pr * sin(pt);
    fprintf('the point ( %d , %d ) in polar is ( %d, %d ) in rect\n\n',...
        pt, pr, recx, recy);
    
    % now to the big boy way
    display('big-boy');
    [pt, pr] = cart2pol(recx, recy);
    [recx, recy] = pol2cart(pt, pr);
    fprintf('\nthe point ( %d , %d ) in rect is ( %d, %d ) in polar\n',...
        recx, recy, pt, pr);
    fprintf('the point ( %d , %d ) in polar is ( %d, %d ) in rect\n\n',...
        pt, pr, recx, recy);
    %{
    figure(1);   % just fucking around
    hold on;
    pax = polaraxes;
    th = 0:.01:2*pi;
    rh = sin(2*th).*cos(2*th);
    polarplot(th, rh);
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 12;
    
    figure(2); % just fucking around
    hold on;
    pax = polaraxes;
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 12;
    polarplot(1:10);
    %}
    figure(3); % the point in polar
    hold on;
    pax = polaraxes;
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 12;
    rlim(pax,[0,pr+1]);
    polarplot(pt, pr, 'r.','markersize',30);
    hold on;
    polarplot([0, pt], [0, pr], 'g:', 'linewidth', 3);
    
    figure(4)
    hold on;
    grid on;
    xlabel('X axis');
    ylabel('Y axis');
    xlim([-abs(recx)-1,1+abs(recx)]);
    ylim([-abs(recy)-1,1+abs(recy)]);
    xax = linspace(-abs(recx)-1, abs(recx)+1, 128);
    yax = linspace(-abs(recy)-1, abs(recy)+1, 128);
    plot(xax, 0*xax, 'k', 'linewidth', 1);
    plot(0*yax, yax, 'k', 'linewidth', 1);
    plot(recx, recy, 'b.','markersize',29);
    hold off;
end
    

%----------------------------------------------------------------------------------------------- #2
if sel == 2
    rx = 1;   % rectangular x component, CHANGE
    ry = 2;   % rectangular y component, CHANGE
    rz = 3;   % rectangular z component, CHANGE
    ct = 0;   % cylindrical theta param...known as fi in book
    cr = 0;   % cyclindrical rho param...known as rho in book
    cz = 0;   % cylindrical z component...equiv to z in rectangular
    
    % the (x, y, z) -> (theta, rho, z)
    cr = sqrt( rx^2 + ry^2 );
    ct = atan(ry / rx);
    cz = rz; % easy
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in cylindrical\n',...
        rx, ry, rx, ct, cr, cz);
    % check (theta, rho, z)  ->  (x, y, z) 
    rx = cr * cos(ct);
    ry = cr * sin(ct);
    rz = cz;
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in cylindrical\n',...
        rx, ry, rx, ct, cr, cz);
    
    fprintf('\n\n  now for the big boy way \n\n');
    [rx, ry, rz] = pol2cart(ct, cr, cz);
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in cylindrical\n',...
        rx, ry, rx, ct, cr, cz);
    [ct, cr, cz] = cart2pol(rx, ry, rz);
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in cylindrical\n',...
        rx, ry, rx, ct, cr, cz);
end

%----------------------------------------------------------------------------------------------- #3
if sel == 3
    rx = 4;   % rectangular x component, CHANGE
    ry = 4;   % rectangular y component, CHANGE
    rz = 4;   % rectangular z component, CHANGE
    sa = 0;   % cylindrical theta param...known as fi in book
    se = 0;   % cyclindrical rho param...known as rho in book
    sr = 0;   % cylindrical z component...equiv to z in rectangular
    
    % the (x, y, z) -> (az, el, rh)
    sa = atan(ry / rx);
    se = atan( sqrt(rx^2 + ry^2) / rz );
    sr = sqrt( rx^2 + ry^2 + rz^2); % easy
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in spherical\n',...
        rx, ry, rx, sa, se, sr);
    
    % check (az, el, rh)  ->  (x, y, z) 
    rx = sr * sin(se) * cos(sa);
    ry = sr * sin(se) * sin(sa);
    rz = sr * cos(se);
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in spherical\n',...
        rx, ry, rx, sa, se, sr);
    
    fprintf('\n\n  now for the big boy way \n\n');
    % this is the important part     book way -> MatLab    need other ifs doing matlab -> book
    if se <= (pi/2)
        se = (pi/2) - se;
    end
    if se > (pi/2)
        se = -(se - (pi/2) );
    end
    if sa > pi
        sa = sa - 2*pi;
    end
    %if sa <= pi   then it is good
    [rx, ry, rz] = sph2cart(sa, se, sr);
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in spherical\n',...
        rx, ry, rx, sa, se, sr);
    [sa, se, sr] = cart2sph(rx, ry, rz);
    fprintf('the point ( %d, %d, %d ) in rectangular is ( %d, %d, %d ) in spherical\n',...
        rx, ry, rx, sa, se, sr);
end






%-------------------------------------------------------------------------------------  RE EDUCATION
%{
    syms rx;                % representation of rectangular x component
    assume(rx, 'real');
    syms ry;                % representation of rectangular y component
    assume(ry, 'real');
    syms pt;                % representation of polar theta component
    assume(pt, 'real');
    syms pr;                % representation of polar rho component
    assume(pr, 'real');
    
    % the transformation formulas              if it is not in the limit of atan, angle will be off
    pt = atan(ry / rx);
    pr = sqrt( rx^2 + ry^2 );
    rx = pr * cos(pt);  % do not define them until you are ready to use them
    ry = pr * sin(pt);
%}
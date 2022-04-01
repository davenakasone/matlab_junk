%{
    ch13.3
    #1    |z + 1 - j5| <= (3/2)   
    #2    0 < |z| < 1
    #3    pi < |z - 4 + j2| < 3pi
    #4    -pi < imag(z) < pi
    #5    abs(arg(z)) < pi/4
    #6    real(1/z) < 1
    #7    real(z) >= -1
    #8    |z + i| >= |z-i|
    #10   f(z) = 5z^2 -12z + 3 + j2  @ 4-j3
    #11   f(z) = 1/(1-z)  @ 1-j
    #12   f(z) = (z-2)/(z+2)  @ 0+j8
    #14   f(z) = real(z^2)/abs(z)    f(0) = 0 if z = 0      it is continuous and imag always = 0
    #15   f(z) = |z|^2 * imag(1/z)   f(0) = 0 if z = 0      it is continuous
    #16   f(z) = imag(z^2) / |z|^2   f(0) = 0 if z = 0      not continuous
    #17   f(z) = real(z) / (1 - |z| )   f(0) = 0 if z = 0      it is continuous
    #18   f(z) = (z - j) / (z + j)  @ 0 + j      derivative
    #19   f(z) = (z - j4)^8  @ 3 + j4            derivative
    #20   f(z) = ((3/2)z + j2) / (j3z - 4)  @ any z ... it is always 0
    #21   f(z) =  j(1-z)^n    @ z = 0 + j0 "0"   power rule a lot easier       !!! matlab takes z straight
    #22   f(z) = (j*z^3 + 3*z^2)^3        @  z = 0 + j2 .... 0  power->chain best
    #23   f(z) = z^3 / (z + j)^3        @  z = 0 + j   power->chain best
%}

clc; 
clf;
close all;
clearvars;
select = 23; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice

%----------------------------------------------------------------------------------------------- 1
if select == 1         % circle on -1 + 5j , r = 1.5
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = abs(z + 1 - j*5);
    plane = zeros(length(x), length(y)) + 3/2; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
       
    surf(x, y, plane, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('|z + 1 - j5| <= (3/2)');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [1.5, 1.5]);
    clabel(a, h);
    hold off;
end


%----------------------------------------------------------------------------------------------- 2
if select == 2         % circle on 0,0 , r = 1
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = abs(z);
    plane1 = zeros(length(x), length(y)) + 1; % must be a mesh
    plane2 = zeros(length(x), length(y)); % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('1 < |z| < 1');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [0, 1]);
    clabel(a, h);
    hold off;
end


%----------------------------------------------------------------------------------------------- 3
if select == 3         % 2 circles...
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = abs(z);
    plane1 = zeros(length(x), length(y)) + pi; % must be a mesh
    plane2 = zeros(length(x), length(y)) + 3*pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('pi < |z - 4 + j2| < 3pi');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [pi, 3*pi], 'r');
    clabel(a, h);
    title('pi < |z - 4 + j2| < 3pi');
    hold off;
end

%----------------------------------------------------------------------------------------------- 4
if select == 4         % any real (x) ...imag(y) is simply [-pi, pi]
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = imag(z);
    plane1 = zeros(length(x), length(y)) - pi; % must be a mesh
    plane2 = zeros(length(x), length(y)) + pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('-pi < imag(z) < pi');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [-pi, pi], 'r');
    clabel(a, h);
    title('-pi < imag(z) < pi');
    hold off;
end


%----------------------------------------------------------------------------------------------- 5
if select == 5         % arctan gives away the y=x, y=-x
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = abs(angle(z));
    plane1 = zeros(length(x), length(y)) + pi/4; % must be a mesh
    %plane2 = zeros(length(x), length(y)) + pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    %surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('abs(arg(z)) < pi/4');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [pi/4, pi/4], 'r');
    clabel(a, h);
    title('abs(arg(z)) < pi/4');
    hold off;
end


%----------------------------------------------------------------------------------------------- 6
if select == 6         % pretty much everything except the singularity
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = real(1./z);
    plane1 = zeros(length(x), length(y)) + 1; % must be a mesh
    %plane2 = zeros(length(x), length(y)) + pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    %surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('real(1/z) < 1');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [1, 1], 'r');
    clabel(a, h);
    title('real(1/z) < 1');
    hold off;
end

%----------------------------------------------------------------------------------------------- 7
if select == 7         % just think about the vertical line it makes...
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = real(z);
    plane1 = zeros(length(x), length(y)) - 1; % must be a mesh
    %plane2 = zeros(length(x), length(y)) + pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    %surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('real(z) >= -1');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, f, [-1, -1], 'r');
    clabel(a, h);
    title('real(z) >= -1');
    hold off;
end


%----------------------------------------------------------------------------------------------- 8
if select == 8         % a few ways to do it, algebra to 0, function comparison, ect  anything y>0
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = abs(z+j*1);
    g = abs(z-j*1);
    plane1 = zeros(length(x), length(y)) +(f-g); % must be a mesh
    %plane2 = zeros(length(x), length(y)) + pi; % must be a mesh
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, f, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, g, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, plane1, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    %surf(x, y, plane2, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    title('|z + i| >= |z-i|');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    [a, h] = contour(x, y, plane1, [0, 0], 'r');
    clabel(a, h);
    title('|z + i| >= |z-i|');
    hold off;
end

%----------------------------------------------------------------------------------------------- 10
if select == 10         % f(z) = 5z^2 -12z + 3 + j2  @ 4-j3
    syms u; % helpers
    px = 4;
    py = -3;
    u0 = (px + j*py);
    g(u) = 5*(u^2)-u*12+3+j*2;
    gmag = abs(eval(g(u0)));
    gi = imag(eval(g(u0)));
    gr = real(eval(g(u0)));
    
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = 5*(z.^2)-z.*12+3+j*2;
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, gmag, 'k.', 'markersize', 20);
    plot3(px, py, gi, 'k.', 'markersize', 20);
    plot3(px, py, gr, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 4-j3');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
    fprintf('real part of f(%d , %d) = %d\n', px, py, gr);
    fprintf('imag part of f(%d , %d) = %d\n', px, py, gi);
end
    

%----------------------------------------------------------------------------------------------- 11
if select == 11         % f(z) = 1/(1-z)  @ 1-j
    syms u; % helpers
    px = 1;  %CHANGE
    py = -1;  %CHANGE
    u0 = (px + j*py);
    g(u) = 1/(1-u); %CHANGE
    gmag = abs(eval(g(u0)));
    gi = imag(eval(g(u0)));
    gr = real(eval(g(u0)));
    
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = (1-z).^-1;  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, gmag, 'k.', 'markersize', 20);
    plot3(px, py, gi, 'k.', 'markersize', 20);
    plot3(px, py, gr, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 1-j'); % CHANGE
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
    fprintf('real part of f(%d , %d) = %d\n', px, py, gr);
    fprintf('imag part of f(%d , %d) = %d\n', px, py, gi); 
end
    

%----------------------------------------------------------------------------------------------- 12
if select == 12         % f(z) = (z-2)/(z+2)  @ 0+j8
    syms u; % helpers
    px = 0;  %CHANGE
    py = 8;  %CHANGE
    u0 = (px + j*py);
    g(u) = (u-2)/(u+2); %CHANGE
    gmag = abs(eval(g(u0)));
    gi = imag(eval(g(u0)));
    gr = real(eval(g(u0)));
    
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = (z-2)./(z+2);  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, gmag, 'k.', 'markersize', 20);
    plot3(px, py, gi, 'k.', 'markersize', 20);
    plot3(px, py, gr, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+8j'); % CHANGE
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-10, 10, 512);
    imagAxis = linspace (-10, 10, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
    fprintf('real part of f(%d , %d) = %d\n', px, py, gr); % 15/17
    fprintf('imag part of f(%d , %d) = %d\n', px, py, gi); % 8/17
end
    

%----------------------------------------------------------------------------------------------- 14
if select == 14         % f(z) = real(z^2)/abs(z)    f(0) = 0 if z = 0   it is continuos
    % the limit of f(z) as z->0 needs to equal 0 for this to be continuous
   
    px = 0;  
    py = 0; 
    u0 = (px + j*py);
    
    syms u; % helpers / feelers
    g(u) = real(u^2)/abs(u);     % CHANGE
    lim0 = limit(g, u, u0);
    lim0r = limit(g, u, u0, 'right');
    lim0l = limit(g, u, u0, 'left');
    fprintf('abs limit = %d \n', lim0);
    fprintf('RH limit = %d \n', lim0r);
    fprintf('LH limit = %d \n', lim0l);
    if (lim0 == 0) & (lim0r == 0) & (lim0l == 0)
        fprintf('f(z) is probably continous because limit = 0 , in most directions\n');
        fprintf('\tcheck the polar just to be safe\n');
    end
    
    
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = real(z.^2)./abs(z);  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, 0, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+0j');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-3, 3, 512);
    imagAxis = linspace (-3, 3, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
end
    

%----------------------------------------------------------------------------------------------- 15
if select == 15         % f(z) = |z|^2 * imag(1/z)   f(0) = 0 if z = 0      it is continuous
    % the limit of f(z) as z->0 needs to equal 0 for this to be continuous
   
    px = 0;  
    py = 0; 
    u0 = (px + j*py);
    
    syms u; % helpers / feelers
    g(u) = ((abs(u))^2)*imag(1/u);     % CHANGE
    lim0 = limit(g, u, u0);
    lim0r = limit(g, u, u0, 'right');
    lim0l = limit(g, u, u0, 'left');
    fprintf('abs limit = %d \n', lim0);
    fprintf('RH limit = %d \n', lim0r);
    fprintf('LH limit = %d \n', lim0l);
    if (lim0 == 0) & (lim0r == 0) & (lim0l == 0)
        fprintf('f(z) is probably continous because limit = 0 , in most directions\n');
        fprintf('\tcheck the polar just to be safe\n');
    end
    
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = ((abs(z)).^2).*imag(1./z);  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, 0, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+0j');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-3, 3, 512);
    imagAxis = linspace (-3, 3, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
end
    

%----------------------------------------------------------------------------------------------- 16
if select == 16         %  #16   f(z) = imag(z^2) / |z|^2   f(0) = 0 if z = 0      not continuous
    % the limit of f(z) as z->0 needs to equal 0 for this to be continuous
   
    px = 0;  
    py = 0; 
    u0 = (px + j*py);
    
    syms u; % helpers / feelers
    g(u) = (imag(u^2))/((abs(u))^2);     % CHANGE
    lim0 = limit(g, u, u0);
    lim0r = limit(g, u, u0, 'right');
    lim0l = limit(g, u, u0, 'left');
    fprintf('abs limit = %d \n', lim0);
    fprintf('RH limit = %d \n', lim0r);
    fprintf('LH limit = %d \n', lim0l);
    if (lim0 == 0) & (lim0r == 0) & (lim0l == 0)
        fprintf('f(z) is probably continous because limit = 0 , in most directions\n');
        fprintf('\tcheck the polar just to be safe\n');
    end
    
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = (imag(z.^2))./((abs(z)).^2);  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, 0, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+0j');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-3, 3, 512);
    imagAxis = linspace (-3, 3, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
end
   

%----------------------------------------------------------------------------------------------- 17
if select == 17         %  #17   f(z) = real(z) / (1 - |z| )   f(0) = 0 if z = 0      it is continuous
    % the limit of f(z) as z->0 needs to equal 0 for this to be continuous
   
    px = 0;  
    py = 0; 
    u0 = (px + j*py);
    
    syms u; % helpers / feelers
    g(u) = real(u) / (1 - abs(u));     % CHANGE
    lim0 = limit(g, u, u0);
    lim0r = limit(g, u, u0, 'right');
    lim0l = limit(g, u, u0, 'left');
    fprintf('abs limit = %d \n', lim0);
    fprintf('RH limit = %d \n', lim0r);
    fprintf('LH limit = %d \n', lim0l);
    if (lim0 == 0) & (lim0r == 0) & (lim0l == 0)
        fprintf('f(z) is probably continous because limit = 0 , in most directions\n');
        fprintf('\tcheck the polar just to be safe\n');
    end
    
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = real(z) ./ (1 - abs(z));  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, 0, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+0j');
    hold off;
    
    figure(2);
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-3, 3, 512);
    imagAxis = linspace (-3, 3, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 7, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 7, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
end


%----------------------------------------------------------------------------------------------- 18
if select == 18         %  #18   f(z) = (z - j) / (z + j)  @ 0 + j
    
    px = 0; % change here
    py = 1; % change here
    
    syms x real;
    syms y real;
    z = (x + j*y);
    f = (z - j*1)/(z + j*1);
    fu = real(f);
    fux = diff(fu, x);
    fuy = diff(fu, y);
    fv = imag(f);
    fvx = diff(fv, x);
    fvy = diff(fv, y);
    
    df_x = (fux + j*fvx); % f'(z) = ux + jvx
    df_y = (fvy - j*fuy); % f'(z) = vy - juy
    fprintf("\nf'(z) = %s\n", df_x);
    fprintf("f'(z) = %s\n", df_x);
    dfxp = subs(df_x, [x, y], [px, py]);
    dfxy = subs(df_y, [x, y], [px, py]);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxp);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxy);   
end

%----------------------------------------------------------------------------------------------- 19
if select == 19         %  #19   f(z) = (z - j4)^8  @ 3 + j4
    
    px = 3; % change here
    py = 4; % change here
    
    syms x real;
    syms y real;
    z = (x + j*y);
    f = (z - j*4)^8;
    fu = real(f);
    fux = diff(fu, x);
    fuy = diff(fu, y);
    fv = imag(f);
    fvx = diff(fv, x);
    fvy = diff(fv, y);
    
    df_x = (fux + j*fvx); % f'(z) = ux + jvx
    df_y = (fvy - j*fuy); % f'(z) = vy - juy
    fprintf("\nf'(z) = %s\n", df_x);
    fprintf("f'(z) = %s\n", df_x);
    dfxp = subs(df_x, [x, y], [px, py]);
    dfxy = subs(df_y, [x, y], [px, py]);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxp);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxy);
end


%----------------------------------------------------------------------------------------------- 20
if select == 20         %  #20   f(z) = ((3/2)z + j2) / (j3z - 4)  @ any z ... it is always 0
    
    px = 3; % change here
    py = 4; % change here
    
    syms x real;
    syms y real;
    z = (x + j*y);
    f = ((3/2)*z+j*2)/(j*3*z-4);
    fu = real(f);
    fux = diff(fu, x);
    fuy = diff(fu, y);
    fv = imag(f);
    fvx = diff(fv, x);
    fvy = diff(fv, y);
    
    df_x = (fux + j*fvx); % f'(z) = ux + jvx
    df_y = (fvy - j*fuy); % f'(z) = vy - juy
    fprintf("\nf'(z) = %s\n", df_x);
    fprintf("f'(z) = %s\n", df_x);
    dfxp = subs(df_x, [x, y], [px, py]);
    dfxy = subs(df_y, [x, y], [px, py]);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxp);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxy);
end


%----------------------------------------------------------------------------------------------- 21
if select == 21         %  #21  f(z) =  j(1-z)^n    @ z = 0 + j0 "0"   power rule a lot easier
    
    px = 0; % change here
    py = 0; % change here
    
    syms x real;
    syms y real;
    syms n real;
    z = (x + j*y);
    f = j*(1-z)^n;
    fu = real(f);
    fux = diff(fu, x);
    fuy = diff(fu, y);
    fv = imag(f);
    fvx = diff(fv, x);
    fvy = diff(fv, y);
    
    df_x = (fux + j*fvx); % f'(z) = ux + jvx
    df_y = (fvy - j*fuy); % f'(z) = vy - juy
    fprintf("\nf'(z) = %s\n", df_x);
    fprintf("f'(z) = %s\n", df_x);
    dfxp = subs(df_x, [x, y], [px, py]);
    dfxy = subs(df_y, [x, y], [px, py]);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxp);
    fprintf(' the derivative at [ %d , %d ] = %s\n', px, py, dfxy);
    
    % check this out
    syms u;
    syms m;
    g = j*(1-u)^m;
    dg = diff(g, u);
    fprintf("f'(z) = %s\n", dg);
    temp = subs(dg, u, px+j*py);
    fprintf("f'(0) = %s\n", temp);
end

%----------------------------------------------------------------------------------------------- 22
if select == 22   %    #22   f(z) = (j*z^3 + 3*z^2)^3        @  z = 0 + j2   power->chain best
    
    px = 0;
    py = 2;
    pz = px + j*py;
    syms z;
    f = (j*z^3 + 3*z^2)^3;                  %     the fast way to handle complex differentiation
    df = diff(f, z, 1);
    result = subs(df, z, pz);
    fprintf('the answer is %s\n', result);
end

%----------------------------------------------------------------------------------------------- 23
if select == 23   %    #23   f(z) = z^3 / (z + j)^3        @  z = 0 + j   power->chain best
    
    px = 0;
    py = 1;
    pz = px + j*py;
    syms z;
    f = (z^3) / (z + j*1)^3;                  %     the fast way to handle complex differentiation
    df = diff(f, z, 1);
    result = subs(df, z, pz);
    fprintf('the answer is %s\n', result);
end







%------------------------------------------------------------------------------------- RE-EDUCATION:
%{
    px = 0;  % CHANGE
    py = 1;  % CHANGE
    u0 = (px + j*py);
    modu = abs(u0);
    
    
    %{
    
    syms u; % helpers / feelers
    g(u) = (u - j*1) / (u + j*1);     % CHANGE
    lim0 = limit(g, u, u0);
    lim0r = limit(g, u, u0, 'right');
    lim0l = limit(g, u, u0, 'left');
    fprintf('abs limit = %d \n', lim0);
    fprintf('RH limit = %d \n', lim0r);
    fprintf('LH limit = %d \n', lim0l);
    %}
    
    x = linspace(-modu-2, modu+2, 100);
    y = linspace(-modu-2, modu+2, 100);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = (z - j*1) ./ (z + j*1);  % CHANGE
    fmag = abs(f);
    fi = imag(f);
    fr = real(f);
    
    figure(1); % the 3D part
    hold on;
    grid on;
    xlabel('real');
    ylabel('imag');
    view(45,45);
    % set x, y, z axis   x is real, y is imag, z is abs(f), imag(f), real(f) or whatever you put
    rax = linspace(-modu-2, modu+2, 512);
    plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 2);
    iax = linspace(-modu-2, modu+2, 512);
    plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 2);
    zax = linspace(-modu*100, modu*100, 512);
    plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 2);
    % place 3 unit vectors i, j, k for refernce
    plot3(1,0,0,'yo', 'markersize', 6, 'linewidth', 2); % i hat
    plot3(0,1,0,'yo', 'markersize', 6, 'linewidth', 2); % j hat
    plot3(0,0,1,'yo', 'markersize', 6, 'linewidth', 2); % k hat
    % graph information in 3D (functions, points, regions, ect)
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    surf(x, y, fmag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
    plot3(px, py, 0, 'k.', 'markersize', 20);
    title('real (red), imag(blue), abs(green)  and values @ 0+j'); % CHANGE
    hold off;
    
    figure(2); % the 2D part, good to see if real(f) orthogonal to imag(f)
    hold on;
    grid on;
    axis equal;
    xlabel('real');
    ylabel('imag');
    realAxis = linspace (-modu-2, modu+2, 1024);
    imagAxis = linspace (-modu-2, modu+2, 1024);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    % plot functions/points here   change contours if needed
    plot(px, py, 'g.', 'markersize', 20);
    [a, h] = contour(x, y, fr, 5, 'r');
    clabel(a, h);
    [b, h] = contour(x, y, fi, 5, 'b');
    clabel(b, h);
    title('looking for ortho');
    hold off;
%}
    

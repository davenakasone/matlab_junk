%{
    *** j is special, and only should be used for (-1)^.5
            don't use it to index of something will steal it  ***
    
    complex differentiation, cauch-riemann for analytic testing
    for f(z) to be analytical, 
        f(z) = u(x,y) + jv(x,y)       real part and imaginary part
            need to see u_x = v_y   and     u_y = -v_x
                or just see if contours of real(f) and imag(f) are orthogonal
    
    f'(z) = u_x + jv_x
    f'(z) = -ju_y + v_y        equate these, and you get the 2 tests
    
    the function should be analytic before you try too much calculus on it
    use tricks to get around poles, undefined areas, ect

    you also have your tried and true:
        f'(z) = lim dz -> 0     [ f(z + dz) - f(z) ] / dz
            or any other common differentiation rule
                to validate   lim real(f) = real(f)   , lim imag(f) = imag(f)   for points


    #1    f(z) = z^2      passes cauchy-reimann       using explicit properties of real/imag part
    #2    f(z) = z^2      uses striaght limit...limits to evaluate, limits for derivative
    #3    f(z) = conj(z)  ...fails     quick script
    #4    f(z) = u(x, y) + jv(x,y)   = e^x cosy + je^x siny    13.4 ex 2
    
%}

clc; 
clf;
close all;
clearvars;
select = 4; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice



%----------------------------------------------------------------------------------------------- 1
%  f(z) = z^2      passes cauchy-reimann test, the function is analytic
if select == 1
    setter = 3;                  % CHANGE HERE to play with z-axis
    px = 1;                      % CHANGE HERE    x/real part of test point
    py = 1;                      % CHANGE HERE    y/imag part of test point
    pz = (px + j*py);
    modu = abs(pz);
    syms x real;
    syms y real;
    f = (x + j*y)^2;  % f(z) where z = (x + jy)        CHANGE HERE
    u = real(f); 
    u_x = diff(u, x);
    u_y = diff(u, y);
    v = imag(f);
    v_x = diff(v, x);
    v_y = diff(v, y);
    f_mag = abs(f); % just for the plane on abs(f)
    f_magx = diff(f_mag, x);
    f_magy = diff(f_mag, y);
    nor = [-f_magx, -f_magy, 1];
    norm = (1 + (f_magx)^2 + (f_magy)^2)^(1/2);
    norma = nor./norm;
    normal1 = subs(norma(1), [x,y], [px, py]); % x comp of unit normal vector @ given point
    normal2 = subs(norma(2), [x,y], [px, py]); % y comp
    normal3 = subs(norma(3), [x,y], [px, py]); % z comp
    fprintf('f(z) = %s  \n\n ', f);              
    test1 = subs(u_x, [x, y], [px, py]) - subs(v_y, [x, y], [px, py]);   % test 1,    u_x - v_y = 0
    test2 = subs(u_y, [x, y], [px, py]) + subs(v_x, [x, y], [px, py]);   % test 2,    u_y + v_x = 0
    fprintf('test1, u_x = %s    v_y = %s   \n', u_x, v_y);
    if test1 == 0
        fprintf('test1 was succesful,   u_x = v_y\n\n');
    else
        fprintf('test1 was a failure, u_x is not v_y      the function is not analytic\n\n');
    end
    fprintf('test2, u_y = %s    v_x = %s   \n', u_y, v_x);
    if test2 == 0
        fprintf('test2 was succesful,   u_y = -v_x\n\n');
    else
        fprintf('test2 was a failure, u_y is not -v_x      the function is not analytic\n\n');
    end
    
    if test1 == 0 && test2 == 0
        % analyze the function  if it passed test                  
        syms fzx;
        syms fzy;
        fzx = (u_x + j*v_x);
        fzy = (v_y - j*u_y);
        xVal = subs(fzx, [x, y], [px, py]);
        yVal = subs(fzy, [x, y], [px, py]);
        fVal = subs(f, [x, y], [px, py]);
        fprintf("at [ %d , %d ] f(z) =  %s \n", px, py, fVal);
        fprintf("f'(z) = %s  =  %s\n" , fzx, xVal); % through d/dx
        fprintf("f'(z) = %s  =  %s\n" , fzy, yVal); % through d/dy

        xVec = linspace(-2-modu, 2+modu, 20);
        yVec = linspace(-2-modu, 2+modu, 20);
        [xVec, yVec] = meshgrid(xVec, yVec);
        z = (xVec + j*yVec);
        g = (z.^2);                                           % CHANGE HERE , use vector operators
        g_mag = abs(g);
        gv = imag(g);
        gu = real(g);
        g_phase = angle(g);

        rax = linspace(-modu-2, modu+2, 128);
        iax = linspace(-modu-2, modu+2, 128);
        zax = linspace(-modu*setter, modu*setter, 128);

        figure('Position', [20, 20, 700, 500]); %   for abs(f), real(f), and imag(f) ,  angle(f) if needed
        hold on;
        grid on;
        xlabel('real');
        ylabel('imag');
        view(160,16);
        % set x, y, z axis   x is real, y is imag, z is abs(f), imag(f), real(f) or whatever you put
        xlim([-modu-2, modu+2]);
        ylim([-modu-2, modu+2]);
        zlim([-modu*setter, modu*setter]);
        plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 2);
        plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 2);
        plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 2);
        % place 3 unit vectors i, j, k for reference
        plot3(1,0,0,'yo', 'markersize', 6, 'linewidth', 2); % i hat
        plot3(0,1,0,'yo', 'markersize', 6, 'linewidth', 2); % j hat
        plot3(0,0,1,'yo', 'markersize', 6, 'linewidth', 2); % k hat
        % graph information in 3D (functions, points, regions, ect)
        surf(xVec, yVec, gu, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
        surf(xVec, yVec, gv, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
        surf(xVec, yVec, g_mag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
        fpz = double(abs(subs(f, [x,y],[px,py])));
        plot3(px, py, fpz, 'y.', 'markersize', 20);         % test point
        tit_str = string(f) + " at [" + string(pz) + "] ...  red = real(f) , blue = imag(f), green = abs(f)";
        title(tit_str, 'fontsize', 20); 

        figure('Position', [730, 20, 700, 500]); %   contours orthogonal
        hold on;
        grid on;
        axis equal;
        xlabel('real');
        ylabel('imag');
        plot(rax, 0*rax, 'k', 'linewidth', 1);
        plot(0*iax, iax, 'k', 'linewidth', 1);
        % plot functions/points here   change contours if needed
        plot(px, py, 'g.', 'markersize', 20);
        [a, h] = contour(xVec, yVec, gu, 5, 'r');
        clabel(a, h);
        [b, h] = contour(xVec, yVec, gv, 5, 'b');
        clabel(b, h);
        title('looking for ortho ,  red = real(f)   vs  blue = imag(f)', 'fontsize', 16);

        figure('Position', [20, 510, 700, 500]);   % diff properties
        hold on;
        grid on;
        set(gca,  'color', [.9, .9, .9]);
        xlabel('real');
        ylabel('imag');
        view(570,33);
        % set x, y, z axis   x is real, y is imag, z is abs(f), imag(f), real(f) or whatever you put
        xlim([-modu-2, modu+2]);
        ylim([-modu-2, modu+2]);
        zlim([-modu*setter, modu*setter]);
        plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 2);
        plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 2);
        plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 2);
        % place 3 unit vectors i, j, k for reference
        plot3(1,0,0,'yo', 'markersize', 6, 'linewidth', 2); % i hat
        plot3(0,1,0,'yo', 'markersize', 6, 'linewidth', 2); % j hat
        plot3(0,0,1,'yo', 'markersize', 6, 'linewidth', 2); % k hat
        %diff prep
        gu_x = zeros(length(xVec), length(yVec));  % du/dx  ... d/dx real(f) 
        for i = 1:length(xVec)
            for k = 1:length(yVec)
                gu_x(i, k) = subs(u_x, [x, y] , [xVec(i,k), yVec(i,k)]);
            end
        end
        gu_y = zeros(length(xVec), length(yVec));    % du/dy  ... d/dy real(f)
        for i = 1:length(xVec)
            for k = 1:length(yVec)
                gu_y(i, k) = subs(u_y, [x, y] , [xVec(i,k), yVec(i,k)]);
            end
        end
        gv_x = zeros(length(xVec), length(yVec));    % dv/dx  ... d/dx imag(f)
        for i = 1:length(xVec)
            for k = 1:length(yVec)
                gv_x(i, k) = subs(v_x, [x, y] , [xVec(i,k),yVec(i,k)]);
            end
        end
        gv_y = zeros(length(xVec), length(yVec));    % dv/dy  ... d/dy imag(f)
        for i = 1:length(xVec)
            for k = 1: length(yVec)
                gv_y(i, k) = subs(v_y, [x, y] , [xVec(i,k),yVec(i,k)]);
            end
        end
        tanP = zeros(length(xVec), length(yVec));
        for i = 1:length(xVec)
            for k = 1:length(yVec)
                tanP(i, k) = (normal1/normal3)*(px - xVec(i,k)) + (normal2/normal3)*(py - yVec(i,k)) + fpz;
            end
        end
        dfx = abs((gu_x + j*gv_x)); % the derivative of f(z) through d/dx   = u_x + jv_x
        dfy = abs((gv_y - j*gu_y)); % the derivative of f(z) through d/dy   = v_y - ju_y
        pzx = double(abs(subs(u_x, [x, y], [px, py]) + j*subs(v_x, [x, y], [px, py])));
        pzy = double(abs((subs(v_y, [x, y], [px, py]) - j*subs(u_y, [x, y] , [px, py]))));
        % graph information in 3D (functions, points, regions, ect)
        surf(xVec, yVec, dfx, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
        mesh(xVec, yVec, dfy, 'edgecolor', 'g');
        surf(xVec, yVec, tanP, 'FaceColor', 'g', 'FaceAlpha', .2, 'EdgeColor', 'none');
        surf(xVec, yVec, g_mag, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none');
        plot3(px, py, pzx, 'y.', 'markersize', 20);         % test point @ abs(dfx)
        plot3(px, py, pzy, 'y*', 'markersize', 10);         % test point @ abs(dfy)
        plot3(px, py, fpz, 'k.', 'markersize', 15);         % test point @ abs(f)
        tit_str =  "red  f'(z) through x    vs    mesh f'(x) through y ";
        title(tit_str, 'fontsize', 20); 
        hold off;
%
    else
        fprintf('do not continue further\n'); % the test has failed at least one of two conditions
    end
end
    


%----------------------------------------------------------------------------------------------- 2
if select == 2         % same f(z) = z^2 , just that using limits
    px = 1;            % CHANGE
    py = 1;            % CHANGE
    u0 = (px + j*py);
   
    syms x real; % helpers / feelers   ... keep symbols real unless going to (0 + j0)  diff() takes real only
    syms y real;
    syms del_x real;
    syms del_y real;
    u = (x + j*y);
    g = (u^2);     % CHANGE
    gu = real(g);
    gux = diff(gu, x);
    guy = diff(gu, y);
    gv = imag(g);
    gvx = diff(gv, x);
    gvy = diff(gv, y);
    
    dgx = (gux + j*gvx);
    dgy = (gvy - j*guy);
    fprintf("cauchy/reimann says f'(z) = %s   =  %s \n", dgx, dgy);
end

%{ 
    
    this would be worthless to do any other way
    
    
    syms M;          subs() doesn't appear to change anything
    syms N;
    F = M*N;
    for i = 1:10
        temp = subs(F, [M, N], [1,2]);
        fprintf('%s   = %d  \n', F, temp);
    end
    %}

%----------------------------------------------------------------------------------------------- 3
if select == 3         % same f(z) = conj(z)    should fail  
    syms x real;
    syms y real;
    z = (x + j*y);
    f = conj(z);          % change function here
    fu = real(f);
    fux = diff(fu, x);
    fuy = diff(fu, y);
    fv = imag(f);
    fvx = diff(fv, x);
    fvy = diff(fv, y);
    
    test1 = fux - fvy; % saying      u_x  =  v_y
    test2 = fuy + fvx; % saying      u_y  = - v_x
    if test1 == 0 && test2 == 0
        fprintf('passed Cauchy/Reimann test, function is analytic, u_x  =  v_y  и  u_y  = - v_x\n');
        fprintf('u_x = %s    ,   v_y  = %s \n', fux, fvy);
        fprintf('u_y = %s    ,   v_x  = %s \n', fuy, fvx);
    else
        fprintf('this function failed, do not analyze it\n');
    end
end


%----------------------------------------------------------------------------------------------- 4
if select == 4         % #4    f(z) = u(x, y) + jv(x,y)   = e^x cosy + je^x siny    13.4 ex 2
    syms x real;
    syms y real;
    f = (exp(x)) * ( (cos(y)) + j*(sin(y)) );
    u = real(f);
    v = imag(f);
    
    ux = diff(u, x, 1);
    uy = diff(u, y, 1);
    vx = diff(v, x, 1);
    vy = diff(v, y, 1);
    
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf('\nux = %s  ,  vy = %s    ux - vy = %s\n', ux, vy, test1);
    fprintf('\nuy = %s  ,  vx = %s    uy + vy = %s\n', ux, vy, test2);
    fprintf('\nboth tests passed , function is analytic\n');
end
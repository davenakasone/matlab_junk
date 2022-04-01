%{
    ch 13.2
    ....first few problems are handled with the root solver
    #28  
    #29
    #30
    #31
%}

clc; 
close all;
clear all;
selector = 31; % CHANGE HERE

%---------------------------------------------------------------------------------------------- 28
if selector == 28
    syms z; 
    f(z) = (z^2)-z*(6-j*2)+17-j*6; % CHANGE ME
    eqn = (z^2)-z*(6-j*2)+17-j*6 == 0; % CHANGE ME
    sols = solve(eqn, z);
    num_sols = length(sols);
    fprintf('there appears to be [ %u ] solutions:\n\n', num_sols);
    for i = 1:length(sols)
        pretty(sols(i));
    end
    for i = 1:num_sols
        fprintf('f( %s ) = %s    ....should be 0\n' , sols(i), f(sols(i)));
    end
    Alim = double(abs(max(sols)));
    
    figure(1); % solutions
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-Alim-2, Alim+2]);
    ylim([-Alim-2, Alim+2]);
    title('solutions');
    realAxis = linspace (-Alim-2, Alim+2, 512);
    imagAxis = linspace (-Alim-2, Alim+2, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'b.', 'markersize', 26); % one point for each solution
        label = 'r_' + string(i) + "( " + string(real(sols(i))) + ' , ' + string(imag(sols(i))) + " )";
        text(real(sols(i))+.5,imag(sols(i))+.5,label,'fontweight','bold','fontsize',14);
    end
    hold off;
    
    x = linspace(-Alim, Alim, 500);
    y = linspace(-Alim, Alim, 500);
    n = linspace(-Alim+1, Alim-1, 7);
    [x, y] = meshgrid(x,y);
    u = (x + j*y);
    w = (u.^2)-u.*(6-j*2)+17-j*6; % CHANGE ME
    wi = imag(w);
    wr = real(w);
    xlabel('real');
    ylabel('imag');
    
    figure(2);
    hold on;
    [a1, h1] = contour(x, y, wi, n, 'b');
    clabel(a1, h1);
    [a2, h2] = contour(x, y, wr, n, 'r');
    clabel(a2, h2);
    title('contours');
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'g.', 'markersize', 26); % one point for each solution
    end
    hold off;
    
    figure(3); % only need abs because it is a solution to make f(z) == 0
    hold on;
    grid on;
    x1 = linspace(-Alim, Alim, 500);
    y1 = linspace(-Alim, Alim, 500);
    [x1, y1] = meshgrid(x1,y1);
    u1 = (x1 + j*y1);
    w1 = abs((u1.^2)-u1.*(6-j*2)+17-j*6); % CHANGE ME
    xlabel('real');
    ylabel('imag');
    surf(x1,y1,w1, 'Edgecolor', 'none');
    colormap('jet');
    colorbar;
    view(60,5);
    title('abs(f(z))');
    points = zeros(length(sols), 3);
    for i = 1:length(sols)
        points(i,1) = double(real(sols(i)));
        points(i,2) = double(imag(sols(i)));
        points(i,3) = eval(f(sols(i)));
        plot3(points(i,1), points(i,2), points(i,3), 'k.', 'markersize', 40);
        label1 = "( " + string(points(i,1)) + " , " + string(points(i,2)) + " , " + string(points(i,3)) + " )";
        text(points(i,1)+.5, points(i,2)+.5, points(i,3)-.5, label1,'fontweight','bold','fontsize',14);
    end
    hold off;
end


%---------------------------------------------------------------------------------------------- 29
if selector == 29
    syms z; 
    f(z) = (z^2)+z+1-j*1; % CHANGE ME
    eqn = (z^2)+z+1-j*1 == 0; % CHANGE ME
    sols = solve(eqn, z);
    num_sols = length(sols);
    fprintf('there appears to be [ %u ] solutions:\n\n', num_sols);
    for i = 1:length(sols)
        pretty(sols(i));
    end
    for i = 1:num_sols
        fprintf('f( %s ) = %s    ....should be 0\n' , sols(i), f(sols(i)));
    end
    Alim = double(abs(max(sols)));
    
    figure(1); % solutions
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-Alim-2, Alim+2]);
    ylim([-Alim-2, Alim+2]);
    title('solutions');
    realAxis = linspace (-Alim-2, Alim+2, 512);
    imagAxis = linspace (-Alim-2, Alim+2, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'b.', 'markersize', 26); % one point for each solution
        label = 'r_' + string(i) + "( " + string(real(sols(i))) + ' , ' + string(imag(sols(i))) + " )";
        text(real(sols(i))+.5,imag(sols(i))+.5,label,'fontweight','bold','fontsize',14);
    end
    hold off;
    
    x = linspace(-Alim, Alim, 500);
    y = linspace(-Alim, Alim, 500);
    n = linspace(-Alim+1, Alim-1, 7);
    [x, y] = meshgrid(x,y);
    u = (x + j*y);
    w = (u.^2)+u+1-j*1; % CHANGE ME
    wi = imag(w);
    wr = real(w);
    xlabel('real');
    ylabel('imag');
    
    figure(2);
    hold on;
    [a1, h1] = contour(x, y, wi, n, 'b');
    clabel(a1, h1);
    [a2, h2] = contour(x, y, wr, n, 'r');
    clabel(a2, h2);
    title('contours');
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'g.', 'markersize', 26); % one point for each solution
    end
    hold off;
    
    figure(3); % only need abs because it is a solution to make f(z) == 0
    hold on;
    grid on;
    x1 = linspace(-Alim, Alim, 500);
    y1 = linspace(-Alim, Alim, 500);
    [x1, y1] = meshgrid(x1,y1);
    u1 = (x1 + j*y1);
    w1 = abs((u1.^2)+u1+1-j*1); % CHANGE ME
    xlabel('real');
    ylabel('imag');
    surf(x1,y1,w1, 'Edgecolor', 'none');
    colormap('jet');
    colorbar;
    view(60,5);
    title('abs(f(z))');
    points = zeros(length(sols), 3);
    for i = 1:length(sols)
        points(i,1) = double(real(sols(i)));
        points(i,2) = double(imag(sols(i)));
        points(i,3) = eval(f(sols(i)));
        plot3(points(i,1), points(i,2), points(i,3), 'k.', 'markersize', 40);
        label1 = "( " + string(points(i,1)) + " , " + string(points(i,2)) + " , " + string(points(i,3)) + " )";
        text(points(i,1)+.5, points(i,2)+.5, points(i,3)-.5, label1,'fontweight','bold','fontsize',14);
    end
    hold off;
end


%---------------------------------------------------------------------------------------------- 30
if selector == 30
    syms z; 
    f(z) = (z^4)+324; % CHANGE ME
    eqn = (z^4)+324 == 0; % CHANGE ME
    sols = solve(eqn, z);
    num_sols = length(sols);
    fprintf('there appears to be [ %u ] solutions:\n\n', num_sols);
    for i = 1:length(sols)
        pretty(sols(i));
    end
    for i = 1:num_sols
        fprintf('f( %s ) = %s    ....should be 0\n' , sols(i), f(sols(i)));
    end
    Alim = double(abs(max(sols)));
    
    figure(1); % solutions
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-Alim-2, Alim+2]);
    ylim([-Alim-2, Alim+2]);
    title('solutions');
    realAxis = linspace (-Alim-2, Alim+2, 512);
    imagAxis = linspace (-Alim-2, Alim+2, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'b.', 'markersize', 26); % one point for each solution
        label = 'r_' + string(i) + "( " + string(real(sols(i))) + ' , ' + string(imag(sols(i))) + " )";
        text(real(sols(i))+.5,imag(sols(i))+.5,label,'fontweight','bold','fontsize',14);
    end
    hold off;
    
    x = linspace(-Alim, Alim, 500);
    y = linspace(-Alim, Alim, 500);
    n = linspace(-Alim+1, Alim-1, 7);
    [x, y] = meshgrid(x,y);
    u = (x + j*y);
    w = (u.^4)+324; % CHANGE ME
    wi = imag(w);
    wr = real(w);
    xlabel('real');
    ylabel('imag');
    
    figure(2);
    hold on;
    [a1, h1] = contour(x, y, wi, n, 'b');
    clabel(a1, h1);
    [a2, h2] = contour(x, y, wr, n, 'r');
    clabel(a2, h2);
    title('contours');
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'g.', 'markersize', 26); % one point for each solution
    end
    hold off;
    
    figure(3); % only need abs because it is a solution to make f(z) == 0
    hold on;
    grid on;
    x1 = linspace(-Alim, Alim, 500);
    y1 = linspace(-Alim, Alim, 500);
    [x1, y1] = meshgrid(x1,y1);
    u1 = (x1 + j*y1);
    w1 = abs((u1.^4)+324); % CHANGE ME
    xlabel('real');
    ylabel('imag');
    surf(x1,y1,w1, 'Edgecolor', 'none');
    colormap('jet');
    colorbar;
    view(60,5);
    title('abs(f(z))');
    points = zeros(length(sols), 3);
    for i = 1:length(sols)
        points(i,1) = double(real(sols(i)));
        points(i,2) = double(imag(sols(i)));
        points(i,3) = eval(f(sols(i)));
        plot3(points(i,1), points(i,2), points(i,3), 'k.', 'markersize', 40);
        label1 = "( " + string(points(i,1)) + " , " + string(points(i,2)) + " , " + string(points(i,3)) + " )";
        text(points(i,1)+.5, points(i,2)+.5, points(i,3)-.5, label1,'fontweight','bold','fontsize',14);
    end
    hold off;
end


%---------------------------------------------------------------------------------------------- 31
if selector == 31
    syms z; 
    f(z) = (z^4)-j*6*(z^2)+16; % CHANGE ME
    eqn = (z^4)-j*6*(z^2)+16 == 0; % CHANGE ME
    sols = solve(eqn, z);
    num_sols = length(sols);
    fprintf('there appears to be [ %u ] solutions:\n\n', num_sols);
    for i = 1:length(sols)
        pretty(sols(i));
    end
    for i = 1:num_sols
        fprintf('f( %s ) = %s    ....should be 0\n' , sols(i), f(sols(i)));
    end
    Alim = double(abs(max(sols)));
    
    figure(1); % solutions
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-Alim-2, Alim+2]);
    ylim([-Alim-2, Alim+2]);
    title('solutions');
    realAxis = linspace (-Alim-2, Alim+2, 512);
    imagAxis = linspace (-Alim-2, Alim+2, 512);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'b.', 'markersize', 26); % one point for each solution
        label = 'r_' + string(i) + "( " + string(real(sols(i))) + ' , ' + string(imag(sols(i))) + " )";
        text(real(sols(i))+.5,imag(sols(i))+.5,label,'fontweight','bold','fontsize',14);
    end
    hold off;
    
    x = linspace(-Alim, Alim, 500);
    y = linspace(-Alim, Alim, 500);
    n = linspace(-Alim+1, Alim-1, 7);
    [x, y] = meshgrid(x,y);
    u = (x + j*y);
    w = (u.^4)-j*6*(u.^2)+16; % CHANGE ME
    wi = imag(w);
    wr = real(w);
    xlabel('real');
    ylabel('imag');
    
    figure(2);
    hold on;
    [a1, h1] = contour(x, y, wi, n, 'b');
    clabel(a1, h1);
    [a2, h2] = contour(x, y, wr, n, 'r');
    clabel(a2, h2);
    title('contours');
    for i = 1:length(sols)
        plot(real(sols(i)), imag(sols(i)), 'g.', 'markersize', 26); % one point for each solution
    end
    hold off;
    
    figure(3); % only need abs because it is a solution to make f(z) == 0
    hold on;
    grid on;
    x1 = linspace(-Alim, Alim, 500);
    y1 = linspace(-Alim, Alim, 500);
    [x1, y1] = meshgrid(x1,y1);
    u1 = (x1 + j*y1);
    w1 = abs((u1.^4)-j*6*(u1.^2)+16); % CHANGE ME
    xlabel('real');
    ylabel('imag');
    surf(x1,y1,w1, 'Edgecolor', 'none');
    colormap('jet');
    colorbar;
    view(60,5);
    title('abs(f(z))');
    points = zeros(length(sols), 3);
    for i = 1:length(sols)
        points(i,1) = double(real(sols(i)));
        points(i,2) = double(imag(sols(i)));
        points(i,3) = eval(f(sols(i)));
        plot3(points(i,1), points(i,2), points(i,3), 'k.', 'markersize', 40);
        label1 = "( " + string(points(i,1)) + " , " + string(points(i,2)) + " , " + string(points(i,3)) + " )";
        text(points(i,1)+.5, points(i,2)+.5, points(i,3)-.5, label1,'fontweight','bold','fontsize',14);
    end
    hold off;
end


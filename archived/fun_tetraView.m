function fun_tetraView(funIn, centReal, centImag)
%{
    given f(z) in terms of globals, z = X + 1j*Y
    graphs inputs on xy plane   (z = 0)
    graphs outputs on xy + zht plane (z = whatever you want)
    connecting lines show map
%}
xp = centReal;       
yp = centImag;
titleStr = sprintf('%s',funIn);
global X; % need X and Y as the symbols to z = x + jy
global Y;

    zht = 2;      % change distance input and output complex planes are split
    buf = 5;      % chage to pick square neighborhood around point
    pts = 7;     % change number of points to consider here, don't go too big
        
    xin = linspace(-(abs(xp))-buf, (abs(xp))+buf, pts);
    yin = linspace(-(abs(yp))-buf, (abs(yp))+buf, pts);
    [Xg, Yg] = meshgrid(xin, yin);
    f = zeros(pts, pts);
    for i = 1:pts
        for k = 1:pts
            f(i,k) = subs(funIn, [X, Y], [Xg(i,k), Yg(i,k)] );
        end
    end
    
    maxPx = 0;
    maxPy = 0;
    for i = 1:pts
        for k = 1:pts
            temp = f(i,k);
            if abs(real(temp)) > maxPx
                maxPx = abs(real(temp));
            end
            if abs(imag(temp)) > maxPy
                maxPy = abs(imag(temp));
            end
        end
    end
            
    
    figure(3);
    hold on;
    grid on;
    view(125,10);
    title( titleStr,'fontsize', 16); % ???
    xlabel('real');
    ylabel('imag'); 
    xlim([-maxPx-2, maxPx+2]);
    ylim([-maxPy-2, maxPy+2]);
    zlim([-2      , zht+2  ]);
    xax = linspace(-maxPx-2  , maxPx+2 , 256);
    yax = linspace(-maxPy-2  , maxPy+2 , 256);
    zax = linspace(-1        , zht+1   , 256);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
    plot3(maxPx-.5, 0      , 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot3(0      , maxPy-.5, 0         ,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    
    for i = 1:pts
        for k = 1:pts
            plot3(Xg(i, k), Yg(i, k), 0, 'r.', 'markersize', 10);
            plot3(real(f(i, k)), imag(f(i, k)), zht, 'b.', 'markersize', 10);
            plot3([Xg(i, k); real(f(i, k))], [Yg(i, k); imag(f(i, k))], [0; zht], 'g');
        end
    end
end


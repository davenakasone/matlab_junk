%{
    complex function side-by-side

    given f(z) = u(x,y) + j v(x,y)           u is real part, v is imaginary part    x и y params real
    3D graph only allows:
        x, y, abs(f)
        x, y, imag(f)
        x, y, real(f)
        x, y, phase(f)
        x, y, and some other 1D measurement base on f(z)

    to see the true output of f(z) , you need 4 dimensions
    (x, y) in -> f(x, y) out     where x represents real part, and y represents imaginary part
    just use one 2D graph and put 

    Also see 'tetraview'  a 5-projection object that shows progression of inputs and outputs
    

    #1 familiar f(z) = z^2        just taking it on an input line y = x    do anywhere
    #2 quick look at f(z) = z^2   @ z = 1 + j
    #3 trying the neighborhood approach around z = 1 + j for f(z) = z^2
    #4 trying the neighborhood ...but it is a circle z = 1 + j for f(z) = z^2
 
%}

clc; 
clf;
close all;
clearvars;
select = 4; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice

%----------------------------------------------------------------------------------------------- 1
if select == 1    %  #1 familiar f(z) = z^2        just taking it on an input line y = x    
    zht = 5;      % change distance input and output complex planes are split
    pts = 21;     % change number of points to consider here, don't go too big
    xp = 2;       % CHANGE, input x coord
    yp = 2;       % CHANGE, input y coord
    zp = (xp + j*yp);
    modu = abs(zp);
    syms x real;
    syms y real;
    z = (x + j*y);
    f = z^2;         % CHANGE, function to evaluate
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
        
        xin = linspace(-1-modu, 1+modu, pts);
        yin = xin;
        zin = zeros(pts);
        
        xout = zeros(pts);
        yout = zeros(pts);
        zout = zeros(pts) + zht;
        
        for i = 1:pts
            tempz = subs(z, [x, y], [xin(i), yin(i)]);
            tempf = subs(f, z, tempz);
            xout(i) = real(tempf);
            yout(i) = imag(tempf);
        end
        xmax = 0;
        ymax = 0;
        for i = 1:pts
            if abs(xout(i)) >= xmax
                xmax = abs(xout(i));
            end
            if abs(yout(i)) >= ymax
                ymax = abs(yout(i));
            end
        end
        if xmax < modu
            xmax = modu;
        end
        if ymax < modu
            ymax = modu;
        end
        
        rax = linspace(-xmax-2, xmax+2, 128);
        iax = linspace(-ymax-2, ymax+2, 128);
        zax = linspace(-1, zht+2, 128);
        
        figure('Position', [20, 20, 700, 500]);
        hold on;
        colormap('jet');
        grid on;
        xlabel('real  x');
        ylabel('imag  y');
        view(160,16);
        xlim([-xmax-2, xmax+2]);
        ylim([-ymax-2, ymax+2]);
        zlim([-1, zht+2]);
        plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 1);
        plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 1);
        plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 1);
        plot3(1,0,0,'yo', 'markersize', 6, 'linewidth', 2); % i hat
        plot3(0,1,0,'yo', 'markersize', 6, 'linewidth', 2); % j hat
        plot3(0,0,1,'yo', 'markersize', 6, 'linewidth', 2); % k hat
        title('inputs (x,y) mapped to outputs f(x,y) ... z = x + jy', 'fontsize', 16);
        zf = subs(f, [x, y], [xp, yp]);
        zfx = real(zf);
        zfy = imag(zf);
        plot3([xp; zfx], [yp; zfy], [0; zht], 'k');
        for i = 1:pts
            plot3(xin(i), yin(i), 0, 'r.', 'markersize', 10);
            ztemp = subs(z, [x, y], [xin(i), yin(i)]);
            tempf = subs(f, z, ztemp);
            tempfx = real(tempf);
            tempfy = imag(tempf);
            plot3([xin(i); tempfx], [yin(i); tempfy], [0; zht], 'g');
            plot3(tempfx, tempfy, zht, 'b.', 'markersize', 10);
        end  
    else
        fprintf('this function failed, do not analyze it\n');
    end
end


%----------------------------------------------------------------------------------------------- 2
if select == 2   % #2 quick look at f(z) = z^2   @ z = 1 + j
    px = 1; % CHANGE , x coord point
    py = 1; % CHANGE , y coord point
    pz = (px + j*py);
    syms z;
    f = z^2;  % CHANGE , f(z)
    fprintf('\nthe function is f(z) = %s\n' , f);
    df = diff(f, z, 1);
    fprintf("the first derivative of f(z) is f'(z) = %s\n", df);
    fp = subs(f, z, pz);
    dfp = subs(df, z, pz);
    fprintf("\tf[ %d, %d ] = %s\n", px, py, fp); 
    fprintf("\tf'[ %d, %d ] = %s\n", px, py, dfp);   
end


%----------------------------------------------------------------------------------------------- 3
if select == 3    %  #3  f(z) = z^2     neighborhood approach    
    zht = 2;      % change distance input and output complex planes are split
    buf = 2;      % chage to pick square neighborhood around point
    pts = 11;     % change number of points to consider here, don't go too big
    xp = 1;       % CHANGE, input x coord
    yp = 1;       % CHANGE, input y coord
    zp = (xp + j*yp);
    modu = abs(zp);
    syms x real;
    syms y real;
    z = (x + j*y);
    f = z^2;         % CHANGE, function to evaluate
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
        
        xin = linspace(xp-buf, xp+buf, pts);
        Xg = zeros(pts, pts);
        for i = 1:pts
            for k = 1:pts
                Xg(i, k) = xin(k);
            end
        end
        yin = linspace(yp-buf, yp+buf, pts);
        Yg = zeros(pts, pts);
        for i = 1:pts
            for k = 1:pts
                Yg(i, k) = yin(k);
            end
        end
        Yg = transpose(Yg);
        zin = zeros(pts, pts);
        
        xout = zeros(pts, pts);
        yout = zeros(pts, pts);
        zout = zeros(pts, pts) + zht;
        
        for i = 1:pts
            for k = 1:pts
                tempz = subs(z, [x, y], [Xg(i, k), Yg(i, k)]);
                tempf = subs(f, z, tempz);
                xout(i, k) = real(tempf);
                yout(i, k) = imag(tempf);
            end
        end
        xmax = 0;
        ymax = 0;
        for i = 1:pts
            for k = 1:pts
                if abs(xout(i, k)) >= xmax
                    xmax = abs(xout(i, k));
                end
                if abs(yout(i, k)) >= ymax
                    ymax = abs(yout(i, k));
                end
            end
        end
        if xmax < modu
            xmax = modu;
        end
        if ymax < modu
            ymax = modu;
        end
        
        rax = linspace(-xmax-2, xmax+2, 128);
        iax = linspace(-ymax-2, ymax+2, 128);
        zax = linspace(-1, zht+2, 128);
        
        figure('Position', [20, 20, 700, 500]);
        hold on;
        %colormap('jet');
        grid on;
        xlabel('real  x');
        ylabel('imag  y');
        view(130,23);
        xlim([-xmax-2, xmax+2]);
        ylim([-ymax-2, ymax+2]);
        zlim([-1, zht+2]);
        plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 1);
        plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 1);
        plot3(rax, 0*rax, 0*rax+zht, 'k', 'linewidth', 1);
        plot3(0*iax, iax, 0*iax+zht, 'k', 'linewidth', 1);
        plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 1);
        plot3(xmax+1,0,0,'yo', 'markersize', 6, 'linewidth', 4); % i direction
        plot3(xmax+1,0,zht,'yo', 'markersize', 6, 'linewidth', 4); % i direction
        plot3(0,ymax+1,0,'yo', 'markersize', 6, 'linewidth', 4); % j direction
        plot3(0,ymax+1,zht,'yo', 'markersize', 6, 'linewidth', 4); % j direction
        plot3(0,0,zht+1,'yo', 'markersize', 6, 'linewidth', 4); % k direction
        title('inputs (x,y) mapped to outputs f(x,y) ... z = x + jy', 'fontsize', 16);
        zf = subs(f, [x, y], [xp, yp]);
        zfx = real(zf);
        zfy = imag(zf);
        plot3(xp, yp, 0, 'kx', 'markersize', 10, 'linewidth', 4);
        plot3(zfx, zfy, zht, 'kx', 'markersize', 10, 'linewidth', 4);
        for i = 1:pts
            for k = 1:pts
                plot3(Xg(i, k), Yg(i, k), 0, 'r.', 'markersize', 10);
                plot3([Xg(i, k); xout(i, k)], [Yg(i, k); yout(i, k)], [0; zht], 'g');
                plot3(xout(i, k), yout(i, k), zht, 'b.', 'markersize', 10);
            end
        end  
    else
        fprintf('this function failed, do not analyze it\n');
    end
end



%----------------------------------------------------------------------------------------------- 4
if select == 4    %  #3  f(z) = z^2     neighborhood approach     by circle
    zht = 2;      % change distance input and output complex planes are split
    pts = 9;     % change number of points to consider here, don't go too big
    rings = 9;    % change for rings
    px = 2;       % CHANGE, input x coord  
    py = 3;       % CHANGE, input y coord
    modu = 2;     % CHANGE, distance to circle neighborhood
    
    syms x real;
    syms y real;
    z = (x + 1j*y);
    f = z^2;         % CHANGE, function to evaluate
        
    xin = linspace(-modu, modu, pts);
    Xg = zeros(rings, pts);
    count = 1/rings;
    
    for i = 1:rings
        temp = xin;
        xin = count * xin * i;
        for k = 1:pts
            Xg(i,k) = xin(k) + px;
        end
        xin = temp;
    end
    
    yit = zeros(rings, pts);
    yib = zeros(rings, pts);
    count = 1/rings;
    for i = 1:rings
        rad = count * modu * i;
        for k = 1:pts   
            yit(i,k) = sqrt((rad^2)-(Xg(i,k)-px)^2) + py;
            yib(i,k) = -sqrt((rad^2)-(Xg(i,k)-px)^2) + py;
        end
    end
    
    
    figure('Position', [820, 20, 500, 500]);
    hold on;
    axis equal;
    xlabel('real  x');
    ylabel('imag  y');
    xlim([-2-abs(px)-modu, 2+abs(px)+modu]);
    ylim([-2-abs(py)-modu, 2+abs(py)+modu]);
    plot(px, py, 'ko', 'linewidth', 5);
    for i = 1:rings
        for k = 1:pts
            plot(Xg(i,k),yit(i,k), 'r.');
            plot(Xg(i,k),yib(i,k), 'r.');
        end
    end
        xout1 = zeros(pts, pts);
        xout2 = zeros(pts, pts);
        yout1 = zeros(pts, pts);
        yout2 = zeros(pts, pts);
        zout = zeros(pts, pts) + zht;
        
        for i = 1:rings
            for k = 1:pts
                tempz1 = subs(z, [x, y], [Xg(i, k), yit(i, k)]);
                tempz2 = subs(z, [x, y], [Xg(i, k), yib(i, k)]);
                tempf1 = subs(f, z, tempz1);
                tempf2 = subs(f, z, tempz2);
                xout1(i, k) = real(tempf1);
                yout1(i, k) = imag(tempf1);
                xout2(i, k) = real(tempf2);
                yout2(i, k) = imag(tempf2);
            end
        end
        xmax = 0;
        xmax1 = 0;
        xmax2 = 0;
        ymax = 0;
        ymax1 = 0;
        ymax2 = 0;
        for i = 1:pts
            for k = 1:pts
                if abs(xout1(i, k)) >= xmax1
                    xmax1 = abs(xout1(i, k));
                end
                if abs(xout2(i, k)) >= xmax2
                    xmax2 = abs(xout2(i, k));
                end
                if abs(yout1(i, k)) >= ymax1
                    ymax1 = abs(yout1(i, k));
                end
                if abs(yout2(i, k)) >= ymax2
                    ymax2 = abs(yout2(i, k));
                end
            end
        end
        if xmax1 < xmax2
            xmax = xmax2;
        else
            xmax = xmax1;
        end
        if xmax < modu
            xmax = modu;
        end
        if ymax1 < ymax2
            ymax = ymax2;
        else
            ymax = ymax1;
        end
        if ymax < modu
            ymax = modu;
        end

    rax = linspace(-xmax-2, xmax+2, 128);
    iax = linspace(-ymax-2, ymax+2, 128);
    zax = linspace(-1, zht+2, 128);

    figure('Position', [20, 20, 700, 500]);
    hold on;
    grid on;
    xlabel('real  x');
    ylabel('imag  y');
    view(130,23);
    xlim([-xmax-2, xmax+2]);
    ylim([-ymax-2, ymax+2]);
    zlim([-1, zht+2]);
    plot3(rax, 0*rax, 0*rax, 'k', 'linewidth', 1);
    plot3(0*iax, iax, 0*iax, 'k', 'linewidth', 1);
    plot3(rax, 0*rax, 0*rax+zht, 'k', 'linewidth', 1);
    plot3(0*iax, iax, 0*iax+zht, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 1);
    plot3(xmax+1,0,0,'yo', 'markersize', 6, 'linewidth', 4); % i direction
    plot3(xmax+1,0,zht,'yo', 'markersize', 6, 'linewidth', 4); % i direction
    plot3(0,ymax+1,0,'yo', 'markersize', 6, 'linewidth', 4); % j direction
    plot3(0,ymax+1,zht,'yo', 'markersize', 6, 'linewidth', 4); % j direction
    plot3(0,0,zht+1,'yo', 'markersize', 6, 'linewidth', 4); % k direction
    title('inputs (x,y) mapped to outputs f(x,y) ... z = x + jy', 'fontsize', 16);
    zf = subs(f, [x, y], [px, py]);
    zfx = real(zf);
    zfy = imag(zf);
    plot3(px, py, 0, 'kx', 'markersize', 10, 'linewidth', 4);
    plot3(zfx, zfy, zht, 'kx', 'markersize', 10, 'linewidth', 4);
    for i = 1:rings
        for k = 1:pts
            plot3(Xg(i, k), yit(i, k), 0, 'r.', 'markersize', 10);
            plot3(Xg(i, k), yib(i, k), 0, 'r.', 'markersize', 10);
            plot3([Xg(i, k); xout1(i, k)], [yit(i, k); yout1(i, k)], [0; zht], 'g');
            plot3([Xg(i, k); xout2(i, k)], [yib(i, k); yout2(i, k)], [0; zht], 'g');
            plot3(xout1(i, k), yout1(i, k), zht, 'b.', 'markersize', 10);
            plot3(xout2(i, k), yout2(i, k), zht, 'b.', 'markersize', 10);
        end
    end 
    
end


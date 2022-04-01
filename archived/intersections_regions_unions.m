%{
    pcolor() is nice if you don't want meshz() with view(2)

    |z-j| < 1   is a disk of r = 1
    -(1 + j) < z < (1 + j) 
    center z = j
    
    inequalities can get very detalied
    |z - 2| <= 3 w/ |z - j -j| <=2 and real(z) < imag(z)
    
    make a big genorous grid if you are not sure
    limitations if unbonded ... see unboundedness in advance
%}

select = 1;

%intersections AND   ...2 sets
if select == 1
    x = linspace(-4, 4, 1000);
    y = linspace(-4, 4, 1000);
    [x,y] = meshgrid(x,y);
    z = (x + j*y);
    w = abs(z-j*1); % prep set
    v = abs(z+1*j); % prep set
    ww = (w<=2); % set formed, matrix in complex plane
    vv = (v<=2); % set formed, matrix in complex plane
    tt = (ww & vv); % intersection, logical AND ...1 if satisfied, 0 not
    meshz(x, y, -tt); % negative helps see 
    colormap('gray');
    view(2); % directly above
    hold on;
    x = linspace(-4, 4, 9); % need custom grid, don't use 'grid'
    y = x;
    [x,y] = meshgrid(x,y);
    plot(x, y, 'r');
    plot(y, x, 'r');
    axis equal;
    title('intersection |z-j| <= 2 и |z+j| <= 2');
    xlabel('real');
    ylabel('imag');
end

%unions OR ...2 sets
if select == 2
    x = linspace(-4, 4, 1000);
    y = linspace(-4, 4, 1000);
    [x,y] = meshgrid(x,y);
    z = (x + j*y);
    w = abs(z-j*1); % prep set
    v = abs(z+1*j); % prep set
    ww = (w<=2); % set formed, matrix in complex plane
    vv = (v<=2); % set formed, matrix in complex plane
    tt = (ww | vv); % intersection, logical OR ...1 if satisfied, 0 not
    meshz(x, y, -tt); % negative helps see 
    colormap('gray');
    view(2); % directly above
    hold on;
    x = linspace(-4, 4, 9); % need custom grid, don't use 'grid'
    y = x;
    [x,y] = meshgrid(x,y);
    plot(x, y, 'r');
    plot(y, x, 'r');
    axis equal;
    title('intersection |z-j| <= 2 or |z+j| <= 2');
    xlabel('real');
    ylabel('imag');
end

%intersections AND   ... 3 of them    for union, change AND to OR
if select == 3
    x = linspace(-4, 4, 1000);
    y = linspace(-4, 4, 1000);
    [x,y] = meshgrid(x,y);
    z = (x + j*y);
    w = abs(z-j*1); % prep set
    v = abs(z+1*j); % prep set
    s = abs(z-1); % 3rd set preped
    ww = (w<=2); % set formed, matrix in complex plane
    vv = (v<=2); % set formed, matrix in complex plane
    ss = (s<=1);
    tt = (ww & vv & ss); % intersection, logical AND ...1 if satisfied, 0 not
    meshz(x, y, -tt); % negative helps see 
    colormap('gray');
    view(2); % directly above
    hold on;
    x = linspace(-4, 4, 9); % need custom grid, don't use 'grid'
    y = x;
    [x,y] = meshgrid(x,y);
    plot(x, y, 'r');
    plot(y, x, 'r');
    axis equal;
    title('intersection |z-j| <= 1 и |z+j| <= 1  и |z -1| <= 1');
    xlabel('real');
    ylabel('imag');
end

%pcolor 
if select == 4
    x = linspace(-4, 4, 1000);
    y = linspace(-4, 4, 1000);
    [x,y] = meshgrid(x,y);
    z = (x + j*y);
    w = abs(z-j*1); % prep set
    v = abs(z+1*j); % prep set
    ww = (w<=2); % set formed, matrix in complex plane
    vv = (v<=2); % set formed, matrix in complex plane
    tt = (ww & vv); % intersection, logical AND ...1 if satisfied, 0 not
    h = pcolor(x, y, -tt); % pcolor better 
    
    set(h, 'EdgeColor','none');
    colormap('gray');
    grid;
    axis equal;
    set(gca, 'layer','top'); % to see grid...over plot
    
    hold on;
    x = linspace(-4, 4, 9); % need custom grid, don't use 'grid'
    y = x;
    [x,y] = meshgrid(x,y);
    plot(x, y, 'r');
    plot(y, x, 'r');
    axis equal;
    title('intersection |z-j| <= 2 и |z+j| <= 2');
    xlabel('real');
    ylabel('imag');
end
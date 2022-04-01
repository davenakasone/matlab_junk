%{
    2.2 ex1
    x [-3, 3]
    y [-3, 3]
    f(z) = e^z
    
    'r' red means real

%}


select = 2;

% real(f) = -1, 0, 1 , imag(f) = -1, 0, 1
if select == 1
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    n = linspace (-1, 0, 3);
    [x, y] = meshgrid(x,y);
    z = (x +j*y);
    f = exp(z);
    fr = real(f);
    fi = imag(f);
    
    figure(1);
    [c, h] = contour(x, y, fi, n, 'b');
    clabel(c,h);
    hold on;
    [d, h] = contour(x, y, fr, n, 'r');
    clabel(d,h);
    xlabel('real');
    ylabel('imag');
    grid;
    title('real-ortho-imag');
    
    figure(2);
    surf(x,y,fi, 'FaceColor','b', 'FaceAlpha', .5, 'Edgecolor', 'none');
    hold on;
    surf(x,y,fr, 'FaceColor','r', 'FaceAlpha', .5, 'Edgecolor', 'none');
    xlabel('real');
    ylabel('imag');
    grid;
    title('real(f) и imag(f)');
end

% real(f) = -1, 0, 1 , imag(f) = -1, 0, 1 ... matlab selects contour
if select == 2
    x = linspace(-3, 3, 100);
    y = linspace(-3, 3, 100);
    %n = linspace (-1, 0, 3);
    [x, y] = meshgrid(x,y);
    z = (x +j*y);
    f = exp(z);
    fr = real(f);
    fi = imag(f);
    
    figure(1);
    [c, h] = contour(x, y, fi, 3, 'b'); % fuck that, choose your own
    clabel(c,h);
    hold on;
    [d, h] = contour(x, y, fr, 3, 'r'); % fuck that, choose your own
    clabel(d,h);
    xlabel('real');
    ylabel('imag');
    grid;
    title('real-ortho-imag');
    
    figure(2);
    surf(x,y,fi, 'FaceColor','b', 'FaceAlpha', .5, 'Edgecolor', 'none');
    hold on;
    surf(x,y,fr, 'FaceColor','r', 'FaceAlpha', .5, 'Edgecolor', 'none');
    xlabel('real');
    ylabel('imag');
    grid;
    title('real(f) и imag(f)');
end


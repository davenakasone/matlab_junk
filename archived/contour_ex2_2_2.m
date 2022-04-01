%{
    2.2 ex2
    x [-3, 3]
    y [-3, 3]
    f(z) = abs(e^z^2)
    
    'r' red means real

    degrenerate at asymptote
%}


select = 1;

% f has value of 2 and 5
if select == 1
    x = linspace(-10, 10, 100);
    y = linspace(-10, 10, 100);
    [x, y] = meshgrid(x,y);
    z = (x +j*y);
    f = abs(exp((z.^2)));
    
    figure(1);
    [c, h] = contour(x, y, f, [2, 2], 'b');  % only where function = 2
    clabel(c,h);
    hold on;
    [d, h] = contour(x, y, f, [5, 5], 'r'); % only where function = 5
    clabel(d,h);
    hold on;
    [e, h] = contour(x, y, f, [1, 1], 'g'); % only where function = 1 ...degenerate
    clabel(e,h);
    xlabel('real');
    ylabel('imag');
    grid;
    title('real-ortho-imag');
    
    figure(2);
    surf(x,y,f, 'FaceColor','g', 'FaceAlpha', .5, 'Edgecolor', 'none');
    xlabel('real');
    ylabel('imag');
    grid;
    title('f(z)');
end

%{ 
real(f) = -1, 0, 1 , imag(f) = -1, 0, 1 ... matlab selects contour
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
%}




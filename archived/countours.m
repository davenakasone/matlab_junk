%{
    contor plots are another useful plotting tool
    they are like a regular map with contour lines
    f(z) = u(x,y) + jv(x,y)
    curves for |f(z)|    no need to worry about view on 3D

    any analytic complex function is orthagonal
    but doesnt have to be if f'(z) = 0
    use that hold on; hold off;

    and put 2 different color maps on same plot (hide axes)
%}

x = linspace(-2.5, 2.5, 100);
y = linspace(-2.5, 2.5, 100);
n = linspace(-2, 2, 5);
[x,y] = meshgrid(x,y);
z = (x + j*y);
w = z.^2;

figure(1);
wi = imag(w);
[c, h] = contour(x, y, wi, n, 'b');
%colormap('jet');  % no need
xlabel('x real');
ylabel('y imag');
clabel(c, h);
hold on;
wr = real(w);
[d, h] = contour(x, y, wr, n, 'r');
%colormap('jet'); % no need 
xlabel('x real');
ylabel('y imag');
title('real(z^2) orthogonal to imag(z^2)');
grid;
clabel(d, h);


figure(2);
surf(x,y,wi, 'FaceColor','b', 'FaceAlpha', .5, 'Edgecolor', 'none');
%colormap('jet');
xlabel('x real');
ylabel('y imag');
title('real(z^2) и imag(z^2)');
hold on;
surf(x, y, wr, 'FaceColor','r', 'FaceAlpha', .5, 'Edgecolor', 'none');
%colormap('gray');



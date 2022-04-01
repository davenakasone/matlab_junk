%{

    ch2 ex 3 и 4
    
    x[-2.5, 2.5]
    y[-2.5, 2.5]
    f(z) = cosz - coshz   ...use Euler
    no contour where abs(f) = 0...
    0's of non const analytic function is nowhere 0
%}

x = linspace(-2.5, 2.5, 100);
y = linspace(-2.5, 2.5, 100);
n = linspace(-5, 5, 11);
[x, y] = meshgrid(x, y);
z = (x + j*y);
f = (cos(z)+cosh(z));
f_mag = abs(f);
f_phase = angle(f);
f_real = real(f);  % this is the focus for the probelm...could use any
f_img = imag(f);

figure(1);
[a, h] = contour(x, y, f_real, n, 'b');
clabel(a, h);
hold on;
[b, h] = contour(x, y, f_img, n, 'r'); % see orthogonal
clabel(b, h);
hold on;
[c, h] = contour(x, y, f_mag, n, 'g'); 
clabel(c, h);
xlabel('real');
ylabel('imag');
grid;
title('real(f) contours');

figure(2);
surf(x, y, f_img, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
hold on;
surf(x, y, f_real, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
hold on;
surf(x, y, f_mag, 'FaceColor', 'g', 'FaceAlpha', .5, 'EdgeColor', 'none');
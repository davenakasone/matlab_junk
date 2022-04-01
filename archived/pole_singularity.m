%{
    pole singularity of order N at point z_0
    there is a Laurent expansion in a deleted neighborhood
    lim z -> z_0  (z-z_0)^N f(z) = c(-N)  ... not 0 or inf
    
    representation is f(z) = 1 / (z - z0)^N   g(z)
    g(z) is analytic and non zero at z_0, g(z_0) = c(-N)
    f(z) will blow up at c(-N) / (z-z0)^N
    

    order of pole tells you how fast it rises to inf

    
%}

select = 2;

%f(z) = cosz / z has a simple pole at origin (N = 1)
if select == 1
x = linspace(-1,1, 200);
y = linspace(-1,1,200);
[x,y] = meshgrid(x,y);
z = (x + j*y);
f = cos(z)./(z);
fmag = abs(f);
mesh(x, y, fmag);
colormap('jet');
grid on;
title('mag of f(z) = cosz / z');
xlabel('real');
ylabel('imag');
end



%f(z) = 1 / z(z-2)^2 has a few poles and order is higher (N=2)
if select == 2
x = linspace(-1,4, 150);
y = linspace(-1,1,150);
[x,y] = meshgrid(x,y);
z = (x + j*y);
f = z.*(z-2).^2;
f = 1./f;
fmag = abs(f);

figure(1);
mesh(x, y, fmag);
axis([-.5, 2.5, -1, 1, 0, 30]); % limits x[], y[], z[]
colormap('jet');
view(10,15);
title('mag of f(z) = 1 / z(z-2)^2');
xlabel('real');
ylabel('imag');

figure(2);
x = linspace(-1, 3, 100);
y = linspace(-1, 1, 100);
n = [.5, 1, 2];     % how to select specific values of function
[x, y] = meshgrid(x, y);
z = (x + j*y);
f = z.*(z-2).^2;
f = 1./f;
fmag = abs(f);
[d, h] = contour(x, y, f1mag, n);
axis equal;
grid;
hold on;
clabel(d,h);
title('contour, mag of f(z) = 1 / z(z-2)^2');
xlabel('real');
ylabel('imag');

end

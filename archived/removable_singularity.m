%{
    some function are analytic everywhere but a singularity

    f(z) = 1 / { (z^2 +9)(z-1)^2   has poles z = +/- 3j, and 1
    g(z) = 1 / sin(piz) has infinite poles
    h(z) = e^(1/z) has an isolated essential singularity at z = 0
    
    "Laurent expansion" can have infinite terms
    if the signularities are not isolated, neighborhood will have inf
    there are also non-isolated essential singulars like 1/sin(pi/z)
    
    use branch cuts and L'Hopital to piecewise
%}

% removable singular points  f(z) = sinz / z
x = linspace(-5, 5, 101);
y = linspace(-1, 1, 51);
[x, y] = meshgrid(x, y);
z = (x + j*y);
f = sin(z+eps)./(z+eps); % eps for error, '.' for matrix
                             % no /by 0
                             % z = x + j*y+eps
                             % z = z + (z==0)*eps  lots of ways
                             % x = linspace(-5, 5, 100)
fi = imag(f);
mesh(x, y, fi);
colormap('jet');
grid on;
colorbar;
title('imag(sinz / z)');
xlabel('real');
ylabel('imag');


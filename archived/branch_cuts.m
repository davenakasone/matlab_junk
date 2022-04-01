%{
    isolated essential singularity at z0 ,
    put in series, get principal part
    
    Picard's Theorem says behavior near signularity makes every
    neighborhood lets function take all possible values, with one exeption,
    an infinite number of times

    f(z) = e^(1/z)  z not 0, have to consider all approaches...that means
    real(f), mag(f), imag(f), phase(f), all need to be looked at

    branch cut or branch point singularity in log(z) or z^B B not integer
    do crazy shit also. log(1) could be 0 or j2pi or j2kpi
    z^(1/2) could be +/- j
    
    pick a path not on branch and assign values
    some functions require more than one cut
    now the function is analytic if cut correctly
    
    all cuts with common point mean there is a branch point
    
    you did it correct if function is discontinous on branch
    discontinuity in real or imag indicate branch cut
    matlab uses principal value, not principal branch
    matlab also uses limits to handle cliffs, even though undfined
    matlab looks at everything as exp(logB)
%}

select = 4;

% f(z) = sqrt(z) ... not z^(1/2) , shows matlab's auto select
if select == 1
    x = linspace(-1, 1, 50);
    y = linspace(-1, 1, 50);
    [X, Y] = meshgrid(x, y); % notice capitals
    z = (X + j*Y);
    f = sqrt(z);
    fi = imag(f);
    mesh(X, Y, fi);
    colormap('jet');
    grid on;
    view(45,30);
    title('imag(sqrt(z))');
    xlabel('real');
    ylabel('imag');
end % discontinuos on y = 0, x <= 0

% f = log(z) with contour
if select == 2
    x = linspace(-5, 5, 1000);
    y = x;
    [x,y] = meshgrid(x,y);
    z = (x + j*y);
    f = log(z);
    fi = imag(f);
    del = .01; % experiment to find
    fi(abs(y) < del & x<=0) = nan; % get rid of matlabs auto-range
    % 'not a number' won't be used in function, considered undf
    [g, h] = contour(x, y, fi);
    clabel(g, h);
    grid;
    title('imag(log(z))');
    xlabel('real');
    ylabel('imag');
end

% f(z) = log((z-1)/(z+1))
if select == 3;
    x = linspace(-5, 5, 100);
    y = x;
    [x, y] = meshgrid(x, y);
    z = (x + j*y);
    f = log( (z - 1) ./ (z + 1) );
    fi = imag(f);
    del = .01;
    fi(abs(y) < del & (-1<=x<=1) )=nan;
    [g, h] = contour(x, y, fi, 6); % looking for 6 levels
    clabel(g, h);
    grid;
    title('contour imag(log(z-1/z+1))');
    xlabel('real');
    ylabel('imag');
end

% branch cuts cause problems in real and imaginary portion
% f(z) = z ^ j >>> e^(jlog(z)) >>> e^(-arg(z)) * [cos(log(abs(z)) +jsin(log(abs(z)) ]
if select == 4
    x = linspace(-5, 5, 1000);
    y = x;
    [x, y] = meshgrid(x,y);
    z = (x +j*y);
    f = z.^j;
    del = .01;
    levels = (0:5); % 0, 1, 2, 3, 4, 5 for countour levles
    fi = imag(f);
    fr = real(f);
    
    figure(1);
    fr(abs(y)<del & -1<=x<=1) = nan; % account for branch cut
    [g, h] = contour(x, y, fr, levels, 'r');
    clabel(g, h);
    hold on;
    fi(abs(y)<del & -1<=x<=1) = nan; % account for branch cut
    [G, H] = contour(x, y, fi, levels, 'b');
    clabel(G, H);
    title('imaginary and real for f(z)');  % see they are orthogonal
    xlabel('real');
    ylabel('imag');
    hold off;
    
    figure(2);
    surf(x, y, fi, 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none');
    hold on;
    surf(x, y, fr, 'FaceColor', 'r', 'FaceAlpha', .5, 'EdgeColor', 'none');
    hold on;
    title('imaginary and real for f(z)');  % see they are orthogonal
    xlabel('real');
    ylabel('imag');
    hold off;
    
    display('fuck you');
end
    
    
    
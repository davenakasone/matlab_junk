%{
    even something simple like w(z) = z^2 is hard
    the value of w has a real and imaginary part....so does input  ... that
    is 4D or two complex planes side by side

    w = f(z) = f(x + jy) is good with u(x,y) and v(x,y)
    w = u + jv
    a lot of people like to plot magnitude on z axis and color the phase
    
    a mesh is one solution, but it is rectangular
    meshgrid() takes x and y

    use ndgrid() instead of meshgrid if you need non-uniform
%}

selector = 6; % CHANGE ME


if selector == 1
x = [0, 1, 2]; % x = 0, x = 1, x = 2
y = [3, 5, 6, 7]; % y = 3 ... y = 7    doesn't need even spacing
[x, y] = meshgrid(x,y); % row vectors changed to matrix
display(x);
display(y);
% make the number
z = (x + j*y);
display(z);
% use z values to evaluate w = z^2
w = z.^2;   % . is very important...square each elm vs entire matrix
display(w);
u = real(w);
display(u);
v = imag(w);
display(v);
    % just the imaginary part "v"
    mesh(x,y,v);
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('Imag(z^2)');
    %test/verify
    test = (1 + j*5); % a point in the mesh
    display(imag(test^2));
    display(v(2,2));
end

if selector == 2
    % changing the default view, axis adjustments, real only
    x = linspace(0, 2, 10); 
    y = linspace(0, 7, 20);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    w = z.^2;
    u = real(w);
    mesh(x,y,u);
    view(45,60); % default view(azm, elv)  in degrees
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('Real(z^2)');
end

if selector == 3
    % use surface...don't try and bunch up points
    x = linspace(0, 2, 10); 
    y = linspace(0, 7, 20);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    w = z.^2;
    u = real(w);
    surf(x,y,u);
    view(45,60); % default view(azm, elv)  in degrees
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('Real(z^2)');
end

if selector == 4
    % place figure(i) before you call plotting
    x = linspace(.5, 2, 50);
    y = linspace(.5, 2, 50);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = z + (z.^-1); % f(z) = z + 1/z
    u = real(f);
    v = imag(f);
    
    figure(1);
    meshz(x,y,u);       %meshz puts a plane under function
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('Real(f)');
    
    figure(2);
    meshz(x,y,v);
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('imag(f)'); 
end

if selector == 5
    x = linspace(-2, 2, 40);
    y = linspace(-8, 8, 80);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = exp(z);
    u = real(f);
    v = imag(f);
    w = abs(f);
    t = angle(f); % angle(f) should get the argument
    
    figure(1);
    surfc(x,y,u); % just a contour plot
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('real(f)');
    
    figure(2);
    surfc(x,y,v); 
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('imag(f)');
    
    figure(3);
    surfc(x,y,w); 
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('abs(f)');
    
    figure(4);
    surfc(x,y,t); 
    colormap('jet');
    colorbar;
    xlabel('x real');
    ylabel('y imag');
    title('Arg(f)');
end

if selector == 6
    x = linspace(-5, 5, 500);
    y = linspace(-5, 5, 500);
    [x, y] = meshgrid(x,y);
    z = (x + j*y);
    f = exp(-(z.^2));
    mag = abs(f);
    u = real(f);
    v = imag(f);
    phase = angle(f);
    
    figure(1);
    mesh(x,y,mag);
    xlabel('x real');
    ylabel('y imag');
    title('mag(f)');
    
    figure(2);
    mesh(x,y,phase);
    xlabel('x real');
    ylabel('y imag');
    title('phase(f)');
    
    figure(3);
    mesh(x,y,u);
    xlabel('x real');
    ylabel('y imag');
    title('real(f)');
    
    figure(4);
    mesh(x,y,v);
    xlabel('x real');
    ylabel('y imag');
    title('imag(f)');
end
    
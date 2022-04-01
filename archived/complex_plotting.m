selector = 3;   % change me


if selector == 1
x = linspace(0,10,1000); % 1000 points between 0 and 10
w = sin(x);  % if w is complex, you will miss all the imaginary stuff
plot(x,w);
end

if selector == 2
z = linspace(-2, 2, 101);
f = z.^(1/3);
plot(z,f);    %Warning: Imaginary parts of complex X and/or Y arguments ignored.
grid;

% because z^(1/3) = |z|^(1/3) < (1/3) prin  + (2/3) 2k pi       k = 0, 1, 2
% matlab only takes principal root
end

if selector == 3
z1 = linspace(-2, 2, 101);
f1 = z1.^(1/3);
plot(z1, real(f1), '-', z1, imag(f1), '*');
grid;
end

if selector == 4
theta = linspace(0,2*pi,100); 
r = 10*exp(j*theta); 
plot(r);               % just takes real x imaginary
% plot(r,'.');   % more variety
axis equal; % nice for circles
end

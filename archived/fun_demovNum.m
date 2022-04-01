function [myRoots] = fun_demovNum(numIn,num, den)
%{
    given a number numIn as z = x + jy                   ...no functions
    given radical num/den    as in    numIn^(num/den)
    roots found and graphed

    solving w^(n/m)
%}

syms k real;

n = num;
m = den;
z = numIn;

v = zeros(1, m+1);      % feeding a row vector of 0's enough to make polynomial w^(m) - z = 0
v(1,1) = 1;             % saying 1 * w^(m)
v(1, m+1) = -z;         % saying last coeff (const) is -z    , now w = z^(1/m)
z_roots = roots(v);     % solve it
myRoots = z_roots;

check_prod = ((-1)^(m-1))*-z;  % compare these
display(check_prod);          % should match 
check_roots = -prod(z_roots);  
display(check_roots);         % should match 
check_sum = sum(z_roots);
display(check_sum);           % should be about 0

z_roots = z_roots.^n; % for last part of radical
display(z_roots);
display(z_roots.^(m/n)); % have to check all

% just some plots
modulus = (abs(z))^(n/m);
Smin = -2-abs(z);
Smax = 2+abs(z);

figure;
hold on;
grid on;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title('z and its m roots');
xlabel('Real');
ylabel('Imaginary');
xlim([Smin, Smax]);
ylim([Smin, Smax]);
realAxis = linspace (Smin, Smax, 256);
imagAxis = linspace (Smin, Smax, 256);
plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
plot(real(z), imag(z), 'b.', 'markersize', 16);
h1 = line([0 real(z)], [0 imag(z)], 'linewidth', 1);
set(h1, 'color', 'b');

lbl = 'xxx';
if imag(numIn) < 0
    temp1 = sprintf('%.3f',real(numIn));
    temp2 = sprintf('%.3f', -1*imag(numIn));
    lbl = sprintf('%s - j %s', temp1, temp2);
else
    temp1 = sprintf('%.3f',real(numIn));
    temp2 = sprintf('%.3f', imag(numIn));
    lbl = sprintf('%s + j %s', temp1, temp2);
end
    
   
text(real(z)+.2,imag(z)+.2, lbl,'fontweight','bold','fontsize',14);

for i = 1:m
    plot(real(z_roots(i)), imag(z_roots(i)), 'r*', 'markersize', 16);
    h = line([0 real(z_roots(i))], [0 imag(z_roots(i))], 'linewidth', 1);
    set(h, 'color', 'r');
    %text(real(Z_set(i))+.2,imag(Z_set(i))+.2,'root','fontweight','bold','fontsize',14);
end

param = 0:pi/100:2*pi;
x_u = modulus*cos(param);
y_u = modulus*sin(param);
plot(x_u, y_u,'g');
axis equal;
hold off;

end

%{
finds roots and checks
all values of z^(1/m) can be found by solving w^m - z = 0
w^m - z   should produce m factors (w-w1)(w-w2)...(w-wm) = 0
(w1)(w2)...(wm) = (-1)^(m-1) * z
coeff w^(m-1) = -(w1 + w2 + ... + wm)
can conclude sum is 0 ...only coeff of highest degree and const

    clc; % nice combo to have
    close all;
    clear all;
    
z = (3 + j*4);  % change here
n = 5;  % change here if you have a radical z^(n/m)  uncommon  n can be < 0
m = 7;  % change here           now you have z^(1/m) prepared  m needs to be > 0


v = zeros(1, m+1); % feeding a row vector of 0's enough to make polynomial w^m - z = 0
v(1,1) = 1; % saying 1* w^m
v(1, m+1) = -z; % saying last coeff (const) is -z    , now w = -z^(n/m)
z_roots = roots(v); % solve it

check_prod = ((-1)^(m-1))*-z;  % compare these
display(check_prod);          % should match 
check_roots = -prod(z_roots);  
display(check_roots);         % should match 
check_sum = sum(z_roots);
display(check_sum);           % should be about 0

z_roots = z_roots.^n; % for last part of radical
display(z_roots);
display(z_roots.^(m/n)); % have to check all

% just some plots
modulus = (abs(z))^(n/m);
Smin = -2-abs(z);
Smax = 2+abs(z);

figure;
hold on;
grid on;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title('z and its m roots');
xlabel('Real');
ylabel('Imaginary');
xlim([Smin, Smax]);
ylim([Smin, Smax]);
realAxis = linspace (Smin, Smax, 256);
imagAxis = linspace (Smin, Smax, 256);
plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
plot(real(z), imag(z), 'b.', 'markersize', 16);
h1 = line([0 real(z)], [0 imag(z)], 'linewidth', 1);
set(h1, 'color', 'b');
text(real(z)+.2,imag(z)+.2,'z','fontweight','bold','fontsize',14);

for i = 1:m
    plot(real(z_roots(i)), imag(z_roots(i)), 'r*', 'markersize', 16);
    h = line([0 real(z_roots(i))], [0 imag(z_roots(i))], 'linewidth', 1);
    set(h, 'color', 'r');
    %text(real(Z_set(i))+.2,imag(Z_set(i))+.2,'root','fontweight','bold','fontsize',14);
end

param = 0:pi/100:2*pi;
x_u = modulus*cos(param);
y_u = modulus*sin(param);
plot(x_u, y_u,'g');
axis equal;
hold off;
%}
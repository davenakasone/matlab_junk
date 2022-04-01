% import java.util.*    java can't handle complex numbers...useless   
close all;
clear all;
clc;
% input here  z^(1/n)

z = 216 + j*0;
n = 3;

%Z_set= Stack();
%modulus = abs(z);
%theta = angle(z);
%temp = (modulus^(1/n))*(cos(theta/n)+j*sin(theta/n));
%Z_set.push(temp);
%temp = Z_set.pop();
%display(temp);

% get first root
Z_set = [];
modulus = abs(z);
theta = angle(z);
temp = (modulus^(1/n))*(cos(theta/n)+j*sin(theta/n));
Z_set(1) = temp;

% get other roots, should be n total

for k = 1:n-1
    t = (theta + 2*pi*k)/n;
    temp = (modulus^(1/n))*(cos(t)+j*sin(t));
    Z_set(k+1) = temp;
    k = k + 1;
end

for i = 1:n
    display(Z_set(i));
end

% start graphing

Smin = -2-abs(z);
Smax = 2+abs(z);

figure;
hold on;
grid on;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title('z1 and z2 in complex plane');
xlabel('Real');
ylabel('Imaginary');
xlim([Smin, Smax]);
ylim([Smin, Smax]);
realAxis = linspace (Smin, Smax, 256);
imagAxis = linspace (Smin, Smax, 256);
plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
plot(real(z), imag(z), 'r.', 'markersize', 16);
h1 = line([0 real(z)], [0 imag(z)], 'linewidth', 1);
set(h1, 'color', 'r');
text(real(z)+.2,imag(z)+.2,'z','fontweight','bold','fontsize',14);

%x_cord = 0;
%y_cord = 0;
%centers = [x_cord, y_cord];
%radii = modulus^(1/n);
%viscircles(centers, radii,'b');

param = 0:pi/100:2*pi;
x_u = modulus^(1/n) * cos(param);
y_u = modulus^(1/n) * sin(param);
circle = plot(x_u, y_u,'b');
axis equal;

for i = 1:n
    plot(real(Z_set(i)), imag(Z_set(i)), 'g*', 'markersize', 16);
    h = line([0 real(Z_set(i))], [0 imag(Z_set(i))], 'linewidth', 1);
    set(h, 'color', 'g');
    text(real(Z_set(i))+.2,imag(Z_set(i))+.2,'root','fontweight','bold','fontsize',14);
end

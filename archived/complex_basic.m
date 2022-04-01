%{
j is will give you base imaginary unit (-1)^.5


%}




% I think this just prepares
close all;
clear all;
clc;

% define complex numbers
z1 = 7 - 4j;
z2 = -2 +3j;

% addition ... subtraction is similar
z = z1+z2;
display(z);

% multiplication ... division is similar
z = z1*z2;
display(z);

% conjugate
z1_con = conj(z1);        %conj() is a built in function
display(z1_con);

% magnitude
z1_mag = abs(z1);      % abs() is built in
display(z1_mag);

% angle ...in radians
z1_ang = angle(z1);       % angle() is good to take it polar...even handles quadrant
display(z1_ang);

% plotting
figure;
hold on;
plot(real(z1), imag(z1), 'b.', 'markersize', 16);    % g* is a green star...any way you want
plot(real(z2), imag(z2), 'r.', 'markersize', 16);
grid on;
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlabel('Real');
ylabel('Imaginary');
xlim([-10, 10]);
ylim([-10, 10]);
title('z1 and z2 in complex plane');

% draw axis
realAxis = linspace (-10, 10, 256);
imagAxis = linspace (-10, 10, 256);
plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);

% draw connections point to orgin
h1 = line([0 real(z1)], [0 imag(z1)], 'linewidth', 1);
h2 = line([0 real(z2)], [0 imag(z2)], 'linewidth', 1);
set(h1, 'color', 'b');
set(h2, 'color', 'r');

% label points
text(real(z1)+.2,imag(z1)+.2,'z1','fontweight','bold','fontsize',14); 
text(real(z2)+.2,imag(z2)+.2,'z2','fontweight','bold','fontsize',14);


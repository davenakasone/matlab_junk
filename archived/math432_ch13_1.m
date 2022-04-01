%{
    13.1 problems
        #1 problem 1
        #2 problem 2
        #3 problem 3
        #4 problem 4
        #5 problems 8-15
        #6 problems 16-20
%}

    clc; % nice combo to have
    close all;
    clear all;
select = 6; % CHANGE ME

% #1 powers of j
if select == 1
    z = (0+j*1);
    n = [-5:5];
    
    for i = 1:11
        Z(1,i) = z^(n(i));
    end
    disp(Z.');
    
    figure(1);
    hold on;
    grid on;
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-2, 2]);
    ylim([-2, 2]);
    title('j^n complex plane');
    realAxis = linspace (-10, 10, 256);
    imagAxis = linspace (-10, 10, 256);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    
    plot(real(Z), imag(Z), 'b.', 'markersize', 16); % only goes one of 4 places
    
    param = 0:pi/100:2*pi;
    x_u = cos(param);
    y_u = sin(param);
    plot(x_u, y_u,'g');
    axis equal;
end

% #2 multiply by j, and you rotate by 90°
if select == 2
    z1 = (1 + j*1);
    z2 = (-1 + j*2);
    z3 = (4 - j*3);
    
    z1r = z1*j;
    z2r = z2*j;
    z3r = z3*j;
    
    before = [z1, z2, z3];
    after = [z1r, z2r, z3r];
    
    figure(1);
    hold on;
    grid on;
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-6, 6]);
    ylim([-6, 6]);
    title('z1, z2, and z3 rotated');
    realAxis = linspace (-10, 10, 256);
    imagAxis = linspace (-10, 10, 256);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    axis equal;
    
    plot(real(z1), imag(z1), 'b.', 'markersize', 16); % blue before *j
    h1 = line([0 real(z1)], [0 imag(z1)], 'linewidth', 1);
    set(h1, 'color', 'b');
    text(real(z1)+.2,imag(z1)+.2,'z1','fontweight','bold','fontsize',14); 
    plot(real(z1r), imag(z1r), 'r.', 'markersize', 16); % red after *j
    h1r = line([0 real(z1r)], [0 imag(z1r)], 'linewidth', 1);
    set(h1r, 'color', 'r');
    text(real(z1r)-.5,imag(z1r)-.2,'z1r','fontweight','bold','fontsize',14);
    
    plot(real(z2), imag(z2), 'b.', 'markersize', 16); % blue before *j
    h2 = line([0 real(z2)], [0 imag(z2)], 'linewidth', 1);
    set(h2, 'color', 'b');
    text(real(z2)+.2,imag(z2)-.2,'z2','fontweight','bold','fontsize',14);
    plot(real(z2r), imag(z2r), 'r.', 'markersize', 16); % red after *j
    h2r = line([0 real(z2r)], [0 imag(z2r)], 'linewidth', 1);
    set(h2r, 'color', 'r');
    text(real(z2r)-.5,imag(z2r)-.2,'z2r','fontweight','bold','fontsize',14);
    
    plot(real(z3), imag(z3), 'b.', 'markersize', 16); % blue before *j
    h3 = line([0 real(z3)], [0 imag(z3)], 'linewidth', 1);
    set(h3, 'color', 'b');
    text(real(z3)+.2,imag(z3)+.2,'z3','fontweight','bold','fontsize',14);
    plot(real(z3r), imag(z3r), 'r.', 'markersize', 16); % red after *j
    h3r = line([0 real(z3r)], [0 imag(z3r)], 'linewidth', 1);
    set(h3r, 'color', 'r');
    text(real(z3r)+.2,imag(z3r)+.2,'z3r','fontweight','bold','fontsize',14);
    
    Angles = [angle(z1), angle(z1r); 
              angle(z2), angle(z2r); 
              angle(z3), angle(z3r)];
    display((180/pi)*Angles);
end
    
% #3 (26-j18)/(6-2j)
if select == 3
    z1 = (26 - j*18);
    z2 = (6-j*2);
    z3 = z1 / z2;
    
    figure(1);
    hold on;
    grid on;
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-3, 30]);
    ylim([-30, 3]);
    title('z1, z2, and z3 rotated');
    realAxis = linspace (-30, 30, 256);
    imagAxis = linspace (-30, 30, 256);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    
    plot(real(z1), imag(z1), 'b.', 'markersize', 16); 
    h1 = line([0 real(z1)], [0 imag(z1)], 'linewidth', 1);
    set(h1, 'color', 'b');
    text(real(z1)+.2,imag(z1)+.2,'z1','fontweight','bold','fontsize',14);
    plot(real(z2), imag(z2), 'b.', 'markersize', 16); 
    h2 = line([0 real(z2)], [0 imag(z2)], 'linewidth', 1);
    set(h2, 'color', 'b');
    text(real(z2)+.2,imag(z2)-.2,'z2','fontweight','bold','fontsize',14);
    plot(real(z3), imag(z3), 'r.', 'markersize', 16); 
    h3 = line([0 real(z3)], [0 imag(z3)], 'linewidth', 1);
    set(h3, 'color', 'r');
    text(real(z3)+.2,imag(z3)+.2,'z3','fontweight','bold','fontsize',14);
end

% #4, conj verification
if select == 4
    z1 = (-11 + j*10);
    z1c = conj(z1);
    z2 = (-1 + j*4);
    z2c = conj(z2);
    
    figure(1);
    hold on;
    grid on;
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    xlabel('Real');
    ylabel('Imaginary');
    xlim([-16, 16]);
    ylim([-16, 16]);
    title('z1, z2, conj');
    realAxis = linspace (-20, 20, 256);
    imagAxis = linspace (-20, 20, 256);
    plot(realAxis, 0*realAxis, 'k', 'linewidth', 1);
    plot(0*imagAxis, imagAxis, 'k', 'linewidth', 1);
    axis equal;
    
    plot(real(z1), imag(z1), 'b.', 'markersize', 16); % blue before *j
    h1 = line([0 real(z1)], [0 imag(z1)], 'linewidth', 1);
    set(h1, 'color', 'b');
    text(real(z1)+.2,imag(z1)+.2,'z1','fontweight','bold','fontsize',14); 
    plot(real(z1c), imag(z1c), 'r.', 'markersize', 16); % red after *j
    h1c = line([0 real(z1c)], [0 imag(z1c)], 'linewidth', 1);
    set(h1c, 'color', 'r');
    text(real(z1c)-.5,imag(z1c)-.2,'z1c','fontweight','bold','fontsize',14);
    
    plot(real(z2), imag(z2), 'b.', 'markersize', 16); % blue before *j
    h2 = line([0 real(z2)], [0 imag(z2)], 'linewidth', 1);
    set(h2, 'color', 'b');
    text(real(z2)+.2,imag(z2)-.2,'z2','fontweight','bold','fontsize',14);
    plot(real(z2c), imag(z2c), 'r.', 'markersize', 16); % red after *j
    h2c = line([0 real(z2c)], [0 imag(z2c)], 'linewidth', 1);
    set(h2c, 'color', 'r');
    text(real(z2c)-.5,imag(z2c)-.2,'z2c','fontweight','bold','fontsize',14);
end
    
% #5 problems 8-15
if select ==5
    z1 = (-2+j*11);
    z2 = (2-j*1);
    
    p_8a = z1*z2;
    display(p_8a);
    p_8b = conj(p_8a);
    display(p_8b);
    
    p_9a = real((z1^2));
    display(p_9a);
    p_9b = (real(z1))^2;
    display(p_9b);
    
    p_10a = real(z2^-2);
    display(p_10a);
    p_10b = real(z2^2);
    display(p_10b);
    
    p_11a = ((z1-z2)^2)/16; % it factors the same
    display(p_11a);
    p_11b = ((z1/4)-(z2/4))^2;
    display(p_11b);
    
    p_12a = z1/z2;
    display(p_12a);
    p_12b = z2/z1;
    display(p_12b);
    
    p_13a = (z1 + z2)*(z1 - z2); % should be same
    display(p_13a);
    p_13b = (z1^2)-(z2^2);
    display(p_13b);
    
    p_14a = (conj(z1))/(conj(z2)); % they are equal..
    display(p_14a);
    p_14b = conj(z1/z2);
    display(p_14b);
    
    p15 = 4*(z1+z2)/(z1-z2);
    display(p15);
end

% #6 problems 16-20        passing to a temp then simplify() is the key
if select == 6
    syms x real;
    syms y real;
    z = (x + j*y);
    
    temp = imag(1/z);
    p16a = simplify(temp);
    fprintf('p16a:\n');
    pretty(p16a);
    temp = imag(1/(z^2));
    p16b = simplify(temp);
    fprintf('p16b:\n');
    pretty(p16b);
    
    temp = real(z^4)-(real(z^2))^2;
    p17 = simplify(temp);
    fprintf('p17:\n');
    pretty(p17);
    
    temp = real(((1+j*1)^16)*(z^2));
    p18 = simplify(temp);
    fprintf('p18:\n');
    pretty(p18);
    
    temp = real(z/conj(z));
    p19a = simplify(temp);
    fprintf('p19a:\n');
    pretty(p19a);
    temp = imag(z/conj(z));
    p19b = simplify(temp);
    fprintf('p19b:\n');
    pretty(p19b);
    
    temp = imag( 1 / ( (conj(z))^2)  );
    p20 = simplify(temp);
    fprintf('p20:\n');
    pretty(p20);
    
end
    
    
    
  


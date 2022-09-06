%{
    chapter 4 : time response

    1  :  ex 4.2, 2nd order
    2  :  ex 4.5, 2nd order, find parameters
    3  :  hw, p13a
    4  :  hw, p16b
    5  :  hw, p18
    6  :  hw, p29
    7  :  hw, p31, solve and plot
    8  :  hw, p37
    9  :  hw, p54
    10 :  hw, p56

%}
close all;
clc;
select = 10;

%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    Hs = 200 / (s^2 + 10*s + 200);
    r_den = roots([1, 10, 200]);
end

%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    syms t;
    Gs = 100 / (s^2 + 15*s +100);
    a = 15;
    b = 100;
    wn = sqrt(b)
    zeta = a / (2 * wn);
    Tp = pi / (wn * sqrt(1-zeta^2)); % 0.4750
    os = 100 * exp((-zeta*pi)/sqrt(1-zeta^2)); % 2.84%
    Ts = 4/ (zeta*wn); % 0.5333
    gt = ilaplace(Gs, s, t);
    t_01 = vpasolve(gt == .01, t);
    t_09 = vpasolve(gt == 0.9, t);
    Tr = double(t_09 - t_01); 
end

%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    syms t;
    Gs = 16 / (s^2 + 3*s + 16);
    a = 3;
    b = 16;
    wn = sqrt(b); % 4
    zeta = a / (2 * wn); % 0.3750
    Tp = pi / (wn * sqrt(1 - zeta^2)); % 0.8472
    os = 100 * exp((-pi*zeta)/sqrt(1-zeta^2)); % 28.0597
    Ts = 4 / (wn*zeta); % 2.6667
    g_t = ilaplace(Gs, s, t)
    steady = double(subs(1+g_t, t, 99999));
    t01 = vpasolve(g_t+1==0.1, t);
    t09 = vpasolve(g_t+1==0.9, t);
    Tr = t09-t01;
end

%------------------------------------------------------------------------------------------
if (select == 4)
    os = 8; % percent
    Tp = 10;
    
    zeta = -log(os/100) / sqrt(pi^2 + (log(os/100))^2); % 0.6266
    wn = pi / (Tp*sqrt(1-zeta^2)); % 0.4031
    b = wn^2; % 0.1625
    a = 2*wn*zeta; % 0.5051
    polz = roots([1, a, b]); % -0.2526 +/- 0.3142i
end

%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    syms T;
    A = [2*s^2 + s, -s; -s, s+1];
    Y = [T;0];
    x = A \ Y;
    tht2 = x(2,1);
    Gs = 1 / (2*s^2 + 2*s + 1);
    
    a = 1;
    b = 1/2;
    wn = sqrt(b); % 0.7071
    zeta = a / (2*wn); %0.7071
    os = 100 * exp((-zeta*pi)/sqrt(1-zeta^2)); % 4.3214
    Ts = 4/(zeta*wn); % 8
    Tp = pi / (wn*sqrt(1-zeta^2)); % 6.2832
end

%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    syms t;
    A = [-5, 0; -1, -2];
    B = [3; 1];
    C = [1, 0];
    x0 = [1; 0];
    D = 0;
    
    P = inv(s.*eye(2)-A);
    pretty(P);
    Xs = P * x0 + P * B .* (1/s);
    pretty(simplify(Xs));
    Ys = C*Xs;
    yt = ilaplace(Ys, s, t)
end

%------------------------------------------------------------------------------------------
if (select == 7)
    syms s;
    syms t;
    A = [-3, 1, 0; 0, -6, 1; 0, 0, -5];
    B = [0; 1; 1];
    C = [0, 1, 1];
    x0 = [0; 0; 0];
    P = inv(s.*eye(3)-A);
    Xs = P * x0 + P * B .* (1/s);
    Ys = C*Xs;
    yt = ilaplace(Ys, s, t);
    pretty(yt);
    
    figure;
    grid on;
    axis padded;
    hold on;
    dots = 100;
    tt = linspace(0, 5, dots);
    yy = double(subs(yt, t, tt));
    plot(tt, yy, 'r','LineWidth', 3);
end

%------------------------------------------------------------------------------------------
if (select == 8)
    syms s;
    syms t;
    Gs = (s+(1/2)) / ((s+2)*(s+5));
    Ps = Gs * (1/s);
    pt = ilaplace(Ps, s, t);
    pretty(pt);
    
    [A, B, C, D] = tf2ss([1, 1/2], [1, 7, 10]);
    x0 = [0; 0];
    P = inv(s.*eye(2)-A);
    Xs = P * x0 + P * B .* (1/s);
    Ys = C*Xs;
    yt = ilaplace(Ys, s, t);
    pretty(yt);
    
    figure;
    grid on;
    axis padded;
    hold on;
    dots = 100;
    tt = linspace(0, 5, dots);
    yy = double(subs(yt, t, tt));
    plot(tt, yy, 'r','LineWidth', 3);
end

%------------------------------------------------------------------------------------------
if (select == 9)
    syms z;
    zz = vpasolve(z^2 * (1+pi^2) - 1 == 0, z);
    zeta = zz(2)
    wn = 1 / zeta
    b = wn^2
    a = 2 * wn * zeta
    M = 1/2
    k = b * M
    os = double(100 * exp((-zeta*pi)/sqrt(1-zeta^2)))
end

%------------------------------------------------------------------------------------------
if (select == 10)
    os = 30 / 100;
    zeta = -log(os) / sqrt(pi^2 + (log(os))^2); % 0.3579
    b = 25;
    wn = sqrt(b);
    D = 2*zeta*wn/25; % 0.1431
end
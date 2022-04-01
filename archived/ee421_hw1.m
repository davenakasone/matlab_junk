%{

    EE 421, homework 1

    1:  p 1.13, R = 1meg
    2:  p1.13, R = 1milli

    EE 421, homework 2

    3:  p 1.18, C = 1 pF
    4:  p1.22, C = 1 uF

    EE 421, homework 4
    
    5: p2.16
    6: p2.24

    7: lab2, 10pF delay
%}
clc;
                            select = 7;  % CHANGE CHANGE CHANGE

fprintf("\nselection #%d\n\n", select);
%------------------------------------------------------------------------------------------ #1
if select == 1
    I_s = -1e-6;
    R = 1e6;
    V_in = 1;
    syms V_out;
    
    kcl = ( ((V_out - V_in)/R) + I_s + ((V_out - 0)/R) ) == 0;
    sol = solve(kcl, V_out);
    fprintf("\nV_out =  %f V\n", sol);    % V_out =  1.000000 V
end


%------------------------------------------------------------------------------------------ #2
if select == 2
    I_s = -1e-6;
    R = 1e-3;
    V_in = 1;
    syms V_out;
    
    kcl = ( ((V_out - V_in)/R) + I_s + ((V_out - 0)/R) ) == 0;
    sol = solve(kcl, V_out);
    fprintf("\nV_out =  %f V\n", sol);    % V_out =  0.500000 V
end


%------------------------------------------------------------------------------------------ #3
if select == 3
    f = 200e6;
    syms w;
    w = 2*pi*f;
    C = 1e-12;
    Zc = 1/(1j*w*C);
    V_in = 3.5;
    R = 1000;
    
    V_out = (V_in * Zc) / (R + Zc);
    out_in_mag = 1/sqrt(1+(2*pi*f*R*C)^2);
    V_out_mag = out_in_mag * V_in;
    V_out_ang = angle(V_out);
    
    
    fprintf("V_out_mag = %0.8f V\n", V_out_mag);
    fprintf("V_out_mag = %0.8f V\n", abs(V_out));
    fprintf("V_out_ang = %f rad\n", angle(V_out));
    fprintf("V_out_ang = %f deg\n", rad2deg(angle(V_out)));
end


%------------------------------------------------------------------------------------------ #4
if select == 4
    R = 1000;
    C = 1e-6;
    f = 159;
    w = 2 * pi * f;
    Zc = 1 / (1j * w * C);
    V_in = 2.5;
    
    V_out = V_in * (Zc / (R + Zc));
    fprintf("V_out mag = %f V\n", abs(V_out)); % 1.768628 V
    fprintf("V_out ang = %f rad\n", angle(V_out)); % -0.784911 rad
    fprintf("V_out ang = %f deg\n", rad2deg(angle(V_out))); % -44.972097 deg
end


%------------------------------------------------------------------------------------------ #5
if select == 5
    length = 200e-6;
    width = 2e-6;
    sqR = 1e3;
    bot_area = length * width;
    C_j0bot = 100e-18;
    C_j0sid = 50e-18;
    sf = width;
    % we use .5 um process...
    
    R_total = sqR * (length/width);  % 1.0000e+05 ohms
    C_j0b = C_j0bot * 4;   % 4.0000e-16 per square
    C_j0s = C_j0sid * (2+2+2+2)* 0.5;  % 2.0000e-16
    Cs = C_j0b + C_j0s;  %  6.0000e-16
    t_d = 0.35 * sqR * Cs * 100^2; % 2.1000e-09 s
end


%------------------------------------------------------------------------------------------ #6
if select == 6
    p210_R = 1e6;
    p210_C = 1e-12;
    
    p224_R = 25e6;
    p224_C = 25e-12;
    
    p210_td = .35 * p210_R * p210_C  %  3.5000e-08   35 ns
    p224_td = .35 * p224_R * p224_C  % 2.1875e-04    .22 ms
    ff = p224_td/p210_td
end


%------------------------------------------------------------------------------------------ #7
if select == 7
    C = 10e-12;
    R = 10e3;
    td = (0.7) * R * C;  % 7.0000e-08  70ns
end

fprintf("\n\n\t\t~ ~ ~ PROGRAM COMPLETE ~ ~ ~\n\n");
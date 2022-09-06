%{
    hw6,
    using the simplified circuit diagram of a 3-phase induction motor
%}
clear all;
close all;
clc;
format compact;
select = 4; % g1, g2, g3, g4

syms s;
V_ll = 480;                   % line to line
V_phi = V_ll / sqrt(3);       % per-phase voltage to motor
X = 0.5;                      % X_1 = X_2, the stator and rotor
X_m = 25;                     % the core
R = 0.25;                     % R_1 = R_2, the stator and rotor
nm_nl = 1798.2;               % no load motor speed, rpm
nm_fl = 1746;                 % full load motor speed, rpm
s_nl = 0.001;                 % no load slip
s_fl = 0.03;                  % full load slip
s_step = 0.001;               % vary "s" by this much
f_e = 60;                     % operates at 60 Hz
w_sync = 2 * pi * f_e / 2;    % in radians per second (it is 4-pole), 188.5
polez = 4;                    % 4-pole machine
P_mech = 883;                 % mechanical loss == converted, no load, a constant

% prepare the data
s_ = s_nl:s_step:s_fl;    % no load --> full load
temp = size(s_);
data_ = temp(2);
Rs = R / s;
Rs_ = subs(Rs, s, s_);
n_sync = nm_nl / 1.001;    % in RPM

% get the motor impedance
z_motor = (R + 1j*X) + f_para(1j*X_m, Rs + 1j*X);
z_motor_ = double(subs(z_motor, s, s_));
z_motor_abs = abs(z_motor);
z_motor_abs_ = double(subs(z_motor_abs, s, s_));

% get the angle and pf
theta = angle(z_motor);
theta_ = double(subs(theta, s, s_));
p_f = cos(theta);
p_f_ = double(subs(p_f, s, s_));

% find the currents
I_1 = V_phi / z_motor;
I_1_ = double(subs(I_1, s, s_));
I_1_abs = abs(I_1);
I_1_abs_ = double(subs(I_1_abs, s, s_));

I_2_abs = I_1_abs * X_m / sqrt(Rs.^2 + (X + X_m).^2);
I_2_abs_ = double(subs(I_2_abs, s, s_));

% find the powers
P_in = real(3 * V_phi * conj(I_1));
P_in_ = double(subs(P_in, s, s_));
P_scl = 3 * R * I_1_abs.^2;
P_scl_ = double(subs(P_scl, s, s_));
P_rcl = 3 * R * I_2_abs.^2;
P_rcl_ = double(subs(P_rcl, s, s_));
P_ag = 3 * I_2_abs.^2 * Rs;
P_ag_ = double(subs(P_ag, s, s_));
P_conv = P_ag - P_rcl;
P_conv_ = double(subs(P_conv, s, s_));
P_out = P_conv - P_mech;
P_out_ = double(subs(P_out, s, s_));

% get the efficiency
eta = P_out / P_in;
eta_ = double(subs(eta, s, s_));

% get the torque induced
tau_ind = P_ag / w_sync;
tau_ind_ = double(subs(tau_ind, s, s_));

% solve the motor speed
w_m = P_conv / tau_ind;
w_m_ = double(subs(w_m, s, s_));
f_r = w_m / (2 * pi);
f_r_ = double(subs(f_r, s, s_));
n_m = n_sync - (120*f_r)/polez;
n_m_ = double(subs(n_m, s, s_));

fptr = fopen("hw6_tab.txt", "w");
fprintf(fptr, "n_m(RPM)  |  abs(z)  | angle(z) | abs(I_1) | pf    |  P_in    |  P_scl  |  abs(I_2)  |  P_rcl  |  P_conv  |  P_out  |  t_ind  |  eff\n");
fprintf(fptr, "--------------------------------------------------------------------------------------------------------------------------------------\n");
for ii = 1:1:data_
    fprintf(fptr, "%7.3f   |", n_m_(ii));
    fprintf(fptr, "%7.3f   |", z_motor_abs_(ii));
    fprintf(fptr, "%4.1f      |", rad2deg(theta_(ii)));
    fprintf(fptr, "%5.2f     |", I_1_abs_(ii));
    fprintf(fptr, "%5.3f  |", p_f_(ii));
    fprintf(fptr, "%9.3f   |", P_in_(ii));
    fprintf(fptr, "%9.3f   |", P_scl_(ii));
    fprintf(fptr, "%9.3f |", I_2_abs_(ii));
    fprintf(fptr, "%9.3f |", P_rcl_(ii));
    fprintf(fptr, "%9.3f |", P_conv_(ii));
    fprintf(fptr, "%9.3f |", P_out_(ii));
    fprintf(fptr, "%9.3f |", tau_ind_(ii));
    fprintf(fptr, "%5.3f \n", eta_(ii));
end
fclose(fptr);

%-------------------------------------------------------------------------------------
if select == 1
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title('Motor Impedance', 'fontsize', 20);
    xlabel('s', 'fontsize', 15);
    ylabel('|z| ohms', 'fontsize', 15);
    for ii = 1:1:data_
        plot(s_(ii), z_motor_abs_(ii), 'r^', 'markersize', 15, 'linewidth', 3);
    end
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 2
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title('Power Factor', 'fontsize', 20);
    xlabel('s', 'fontsize', 15);
    ylabel('pf', 'fontsize', 15);
    for ii = 1:1:data_
        plot(s_(ii), p_f_(ii), 'b^', 'markersize', 15, 'linewidth', 3);
    end
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 3
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title('em torque', 'fontsize', 20);
    xlabel('s', 'fontsize', 15);
    ylabel('em torque', 'fontsize', 15);
    for ii = 1:1:data_
        plot(s_(ii), tau_ind_(ii), 'g^', 'markersize', 15, 'linewidth', 3);
    end
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 4
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title('efficiency', 'fontsize', 20);
    xlabel('s', 'fontsize', 15);
    ylabel('eta', 'fontsize', 15);
    for ii = 1:1:data_
        plot(s_(ii), eta_(ii), 'c^', 'markersize', 15, 'linewidth', 3);
    end
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END>  

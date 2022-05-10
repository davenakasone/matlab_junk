%{
    1  :  last year, induction motor problem
    2  :  last year, transmission line
    3  :  last year, transformer
    4  :  last year, generator
%}
format compact;
close all;
clear;
clc;
select = 4;


%-------------------------------------------------------------------------------------
if select == 1
   % 3-phase induction motor
   V_in = 460;
   V_in_ll = sqrt(3) * V_in;
   poles_n = 4;
   w = 188.5; % if 2 pole, use 377 rad/sec
   rpm = 1800; %...4pole is 1800rpm, 2 pole is 3600
   horsez = 20;
   wattz = f_hp2watts(horsez); % convert the power rating
   R_1 = 0.5;

   % below rating, motor runs at 1728 rpm, drawing 12 kW and 5 KVAR
   % total mechanical losses are 300 W   "g_" as given
   g_P = 12e3;
   g_Q = 5e3;
   g_P_loss_mech = 300;
   g_rpm = 1728;

   % 1) find current motor draws
   S_abs = sqrt(g_P^2 + g_Q^2);
   I_draw = S_abs / V_in_ll;
   fprintf("1)  the motor current is:  %0.1f A\n", I_draw);

   % 2) find the airgap power loss
   P_loss_stator = 3 * I_draw^2 * R_1;
   P_loss_ag = g_P - P_loss_stator;
   fprintf("2) P_ag =  %0.1f W\n", P_loss_ag);

   % 3)  total rotor copper loss, P_ag == P_RCL/s 
   n_s = 120 * rpm / poles_n;
   s = (rpm - g_rpm) / rpm;
   P_rcl = P_loss_ag * s;
   fprintf("3)  P_rcl=  %0.1f W\n", P_rcl);

   % 4)  motor output HP, 1HP==746W
   P_out = P_loss_ag - P_rcl - g_P_loss_mech;
   hp_out = P_out / 746;
   fprintf("4)  motor outputs:  %0.1f HP\n", hp_out);

   % 5) efficiency
   eff_eta = P_out / g_P;
   fprintf("5)  %0.1f %%\n", eff_eta*100);
end


%-------------------------------------------------------------------------------------
if select == 2
    % balanaced 3-phase transmission line
    f = 60;
    w = 2 * pi * f;
    V_s_ll = 230e3;
    Z = 5 + 1j*50;
    Y = 1j*0.0008;
    %V_s = V_s_ll / sqrt(3);
    V_s = 130e3; % given
    g_Vr = 124e3; 
    g_theta = deg2rad(-10);
    V_r = g_Vr * exp(1j*g_theta);
    [A, B, C, D] = f_line_med_ABCD(Z, Y);
    f_mdri("A", A, 1);
    f_mdri("B", B, 1);
    f_mdri("C", C, 1);
    f_mdri("D", D, 1);
    fprintf("\n");

    % 1)  find magnitude and phase of current at receiving end
    I_r = (V_s - A*V_r) / B;
    fprintf("1) ");  f_mdri("Ir", I_r, 1);

    % 2) find magintude and phase of current at sending end
    I_s = C*V_r + D*I_r;
    fprintf("2) ");  f_mdri("Is", I_s, 1);

    % 3) the P and Q, receving
    %theta = angle(V_r) - angle(I_r);
    %P_r = abs(V_r) * abs(I_r) * cos(theta);
    %check_Pr = real(V_r * conj(I_r)) - P_r
    S_r = V_r * conj(I_r);
    fprintf("3) "); f_mdri("Sr", S_r, 1/1e6);

    % 4) the P and Q, sending
    S_s = V_s * conj(I_s);
    fprintf("4) "); f_mdri("Ss", S_s, 1/1e6); % it generates reactive power...

    % 5) efficiency
    eta_eff = real(S_r) / real(S_s);
    fprintf("5)  eff=  %0.1f %%\n", 100*eta_eff);

    % 6) if shorted, Vr==
    V_rr = V_s / A;
    fprintf("6) "); f_mdri("Vr", V_rr, 1/1e3);
end


%-------------------------------------------------------------------------------------
if select == 3
    % single phase, 2-winding  given a typical load, secondary rated
    S_rate = 50e3;
    f = 60;
    w = 2 * i * f;
    V_p = 2400;
    V_s = 240;
    a = V_p / V_s;
    z_series_pu = 0.0217 + 1j*0.0434;
    g_P_in = 24.5e3;
    g_P_load = 24e3;
    g_pf_load = 1;

    % 1)  find series impedance, primary side
    z_base = V_p^2 / S_rate;
    z_series = z_base * z_series_pu;
    fprintf("1) "); f_mdri("z_series", z_series, 1);

    % 2) find primary side voltage
    I_s = g_P_load / V_s; % the pf==1
    V_prim = V_p + (I_s/a)*(z_series);
    fprintf("2) "); f_mdri("Vp", V_prim, 1);

    % 3)  find the efficiency
    eff_eta = g_P_load / g_P_in;
    fprintf("3)  eff=  %0.1f %%\n", 100*eff_eta);

    % 4)  find core loss
    P_total = g_P_in - g_P_load;
    P_copper = (I_s/a)^2 * real(z_series);
    P_core = P_total - P_copper;
    fprintf("4) P_core=  %0.1f W\n", P_core);

    % 5) core loss, pu
    R_core =  abs(V_prim)^2 / P_core;
    R_core_pu = R_core / z_base;
    fprintf("5) core loss resistance, pu=  %0.1f\n", R_core_pu);
end


%-------------------------------------------------------------------------------------
if select == 4
    % 3-phase, generator, Y connected
    X_s_pu = 0.2;
    V_phase = 12.47e3;
    S_rate = 100e6;
    pf = 0.8; % lag
    f = 60;
    w = 2 * pi * f;
    theta = -1 * acos(pf);
    % ignore loss

    % 1) internal impedance
    z_base = V_phase^2 / S_rate;
    z_internal = z_base * 1j*X_s_pu;
    fprintf("1) "); f_mdri("z_internal", z_internal, 1);

    % 2) current supplied
    I_A = (S_rate * pf) / (sqrt(3) * V_phase * pf);
    fprintf("2)  generator supplies:  %0.2f kA\n", I_A/1e3);
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END> 
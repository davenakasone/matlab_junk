%{
    field == rotor
    armeture == stator
    induction motor problems

    1  :  p7
    2  : trying the quiz
%}
clear all;
close all;
clc;
format compact;
select = 4;

%-------------------------------------------------------------------------------------
if select == 1
    % Y connected motor, 15 HP rate
    polez = 4;
    f_e = 60;
    w_sync = 2*pi*f_e/2;
    n_sync = 1800;   % you should know this...3600 for p=2, w=377
    V_ll = 208;
    V_phi = V_ll / sqrt(3);
    R_1 = 0.22;
    R_2 = 0.127;
    X_m = 15;
    X_1 = 0.43;
    X_2 = 0.43;
    P_mech = 300;
    P_core = 200;
    s = 0.05;
   
    z_f = f_para( (R_2/s) + 1j*X_2, 1j*X_m);
    f_mdri("z_f", z_f, 1);
    
    I_L = V_phi / ((R_1 + 1j*X_1) + z_f);
    f_mdri("I_L", I_L, 1);
    
    P_scl = 3 * abs(I_L)^2 * R_1;
    f_mdri("P_scl", P_scl, 1);
    
    P_ag = 3 * abs(I_L)^2 * real(z_f);
    f_mdri("P_ag", P_ag, 1);
    
    P_conv = (1-s)*P_ag;
    f_mdri("P_conv", P_conv, 1);
    
    tau_ind = P_ag/w_sync;
    f_mdri("tau_ind", tau_ind, 1);
    
    P_out = P_conv - P_mech - P_core;
    f_mdri("P_out", P_out, 1);
    
    n_m = (1-s)*n_sync
    w_m = n_m * 2 * pi / 60;
    tau_load = P_out / w_m
    
    z_motor = f_z_motor(R_1, X_1, X_m, R_2, X_2, s);
    P_in = 3 * abs(I_L) * V_phi * cos(angle(z_motor));
    
    eta = P_out / P_in
end
  

%-------------------------------------------------------------------------------------
if select == 2
    % 
    polez = 4;
    n_sync = 3600/2; % or 3600, p=2
    w_sync = 377/2; % or 377, p=2
    f_e = 60;
    s = 0.1;
    V_ll = 460;
    %V_phi = V_ll / sqrt(3);
    V_phi = 265.6;
    R_1 = 0.3;
    X_1 = 1;
    R_2 = 0.4;
    X_2 = 1;
    X_m = 30;
    
    P_core = 0; % maybe
    P_mech = 0; % get when off..s=1
    
    z_motor = f_z_motor(R_1, X_1, X_m, R_2, X_2, s);
    f_mdri("z_motor", z_motor, 1); % has the theta...
    fprintf("theta:  %0.3f  ,  pf:  %0.3f\n", rad2deg(angle(z_motor)), cos(angle(z_motor)));
    
    I_1 = conj(V_phi / z_motor);
    f_mdri("I_1", I_1, 1);
    I_1_abs = abs(I_1);
    
    I_2_abs = I_1_abs * X_m / sqrt( (R_2/s)^2 + (X_m + X_2)^2);
    
    P_in = real(3 * V_phi * conj(I_1));
    f_mdri("P_in", P_in, 1);
    
    P_scl = 3 * I_1_abs^2 * R_1;
    f_mdri("P_scl", P_scl, 1);
    
    P_rcl = 3 * I_2_abs^2 * R_2;
    f_mdri("P_rcl", P_rcl, 1);
    
    P_ag =  3 * I_2_abs^2 * (R_2/s);
    f_mdri("P_ag", P_ag, 1);
    
    P_conv = P_ag - P_rcl;
    f_mdri("P_conv", P_conv, 1);
    P_convv = 3 * I_2_abs^2 * R_2 * ((1-s)/s);
    f_mdri("P_convv", P_convv, 1);
    
    P_out = P_conv - P_mech;
    f_mdri("P_out", P_out, 1);
    
    eta = P_out / P_in
    tau_ind = P_ag / w_sync
    w_m = P_conv / tau_ind
    n_m = (1-s)*n_sync
    
end


%-------------------------------------------------------------------------------------
if select == 3
    polez = 4;
    n_sync = 3600/2; % or 3600, p=2
    w_sync = 377/2; % or 377, p=2
    f_e = 60;
    s = 1/30;
    V_ll = 460;
    %V_phi = V_ll / sqrt(3);
    V_phi = 265.6;
    R_1 = 0.3;
    X_1 = 1;
    R_2 = 9;
    X_2 = 1;
    X_m = 30;
    nm_fl = 1740;
    P_rot = 500;
    
    %z_motor = R_1 + 1j*X_m + f_para(1j*X_m, (R_2s + 1j*X_2));
    z_motor = f_z_motor(R_1, X_1, X_m, R_2, X_2, s)
    f_mdri("z_motor", z_motor, 1);
    
    I_1 = V_phi / z_motor;
    f_mdri("I_1", I_1, 1);
    
    I_2 = abs(I_1) * X_m / sqrt(R_2s^2 + (X_m + X_2)^2)
    
    
    
end


%-------------------------------------------------------------------------------------
if select == 4
    
end

%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END>  
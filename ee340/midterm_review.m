%{
    1  : s19, p1 a transformer
    2  : s19, p2 a magnetic circuit
    3  : s21, mag circuit
    4  : s21, trans line
    5  : s21, transformer
    
%}
close all;
clear;
clc;
select = 3;


%-------------------------------------------------------------------------------------
if select == 1
    S_rate = 5e3;
    f = 60;
    w = 2 * pi * f;
    Vp_rate = 480;
    Vs_rate = 120;
    a = Vp_rate / Vs_rate;
    R_unit = 0.05;
    X_unit = 0.1;
    P_loss = 50; % with no load
    i_ex = 312.5e-3;
    
    % series impedance on Primary side...per units
    z_base = Vp_rate^2 / S_rate;
    R_series = z_base * R_unit;
    X_series = z_base * X_unit;
    z_series = R_series + 1j * X_series;
    f_mdri("z_series", z_series, 1);
    
    % find shunt impedance
    S_shunt = Vp_rate * i_ex;
    Q_shunt = sqrt(S_shunt^2 - P_loss^2);
    R_shunt = Vp_rate^2 / P_loss;
    X_shunt = Vp_rate^2 / Q_shunt;
    z_shunt = R_shunt + 1j * X_shunt;
    f_mdri("z_shunt", z_shunt, 1);
    
    % supplied 480V at source, feeds pure resistive load, 3 ohm
    v_source = 480;
    z_load = 3;
    a_vs = v_source * (a^2 * z_load / (a^2 * z_load + z_series));
    v_s = a_vs / a;
    f_mdri("V_s", v_s, 1);
    
    % now get efficiency under above
    i_s = v_s / z_load;
    S_out = v_s * conj(i_s);
    P_out = real(S_out);
    P_copper = real(i_s/a)^2 * real(z_series);
    eta_eff = P_out / (P_out + P_loss + P_copper);
    fprintf("eff=  %0.1f %%\n", 100*eta_eff);
    
    % find voltage regulation
    vr = (Vs_rate - abs(v_s)) / abs(v_s);
    fprintf("VR=  %0.2f %%\n", 100 * real(vr));
    
    % maximize eff, P_core = P_copper, find current, then load
    Is_a = sqrt(P_loss / real(z_series));
    z_eff = Vs_rate / (a*Is_a)
end


%-------------------------------------------------------------------------------------
if select == 2
 
end


%-------------------------------------------------------------------------------------
if select == 3
    % assume linear, no core loss
    f = 60;
    w = 2 * pi * f;
    n_1 = 200;
    n_2 = 20;
    a = n_1 / n_2;
    l_core = 119/100;
    l_air = 1/100;
    a_m = 25e-4;
    mu = 1.2e-3;
    V_rms = 120;
    V_peak = V_rms * sqrt(2);
    
    % flux generated, in RMS
    reluct = (l_core/(a_m*mu)) + (l_air/(a_m*f_perm_mu(1)));
    phi_rms = V_rms / (w * n_1);
    fprintf("\nrms flux made, coil 1:  %0.4f  mWb\n", phi_rms*1e3);
    
    % find RMS value of source current
    B_core = phi_rms / a_m;
    H_core = B_core / mu;
    I_rms = (phi_rms/n_1) * reluct;
    fprintf("I_rms, from source:  %0.1f A\n", I_rms);
    
    % find inductance of energized coil
    L = n_1^2 / reluct; % use this every time....
    fprintf("L of energized coil:  %0.3f mH\n", L*1e3);
    
    % find the rms voltage in the top coil
    v_ind = phi_rms * w * n_2;
    fprintf("top coil gets volts induced:  %0.1f V\n", v_ind);
    
    % recalculate rms current if airgap is gone
    l_path = 120/100; % path is all core
    reluct = (l_core/(a_m*mu));
    phi_rms = V_rms / (w * n_1);
    I_rms = (phi_rms/n_1) * reluct;
    fprintf("I_rms, from source:  %0.1f A\n", I_rms);
    
    % and a 0.5 resistor is added on top coil
    r_load = 0.5;
    E = I_rms * r_load; % ??? wtf
end


%-------------------------------------------------------------------------------------
if select == 4
    % balanced 3-phase line
    f = 60;
    w = 2 * pi * f;
    z_series = 5 + 1j * 50;
    y_shunt = 1j * 0.0008;
    [A, B, C, D] = f_line_med_ABCD(z_series, y_shunt);
    f_mdri("A", A, 1);
    f_mdri("B", B, 1);
    f_mdri("C", C, 1);
    f_mdri("D", D, 1);
    % best to use his given ABCD ?
    A = 0.98 * exp(1j*deg2rad(0));
    B = 50.25 * exp(1j*deg2rad(83.4));
    C = 0.00079 * exp(1j*deg2rad(90));
    D = 0.98 * exp(1j*deg2rad(0));
    V_s = 130e3; % given as phase, not line
    
    % under some load condition:
    V_r = 124e3 * exp(1j*deg2rad(-10));
    % find magnitude and angle of I_r, use V_s = A V_r + B I_r
    I_r = (V_s - A * V_r) / B;
    fprintf("\n");
    f_mdri("I_r", I_r, 1);
    
    % find the per-phase P и Q, delivered
    S_r = V_r * conj(I_r);
    f_mdri("S_r", S_r, 1/1e6);
    
    % mag and phase of I_s
    V_s2 = A * V_r + B * I_r;
    I_s2 = C * V_r + D * I_r;
    f_mdri("I_s", I_s2, 1);
    
    % per-phase P и Q at sending end
    S_s = V_s2 * conj(I_s2);
    f_mdri("S_s", S_s, 1/1e6);
    
    % efficiency:
    P_out = real(S_r);
    P_in = P_out + abs(I_s2)^2 * real(z_series);
    eta_eff = P_out / P_in;
    fprintf("efficiency:  %0.1f %%, making more reactive than consuming\n", 100*eta_eff);
    
    % find charging current, I_s under no load:
    I_s_nl = V_s * C;
    f_mdri("I_s_nl", I_s_nl, 1);
end


%-------------------------------------------------------------------------------------
if select == 5
    % 2-winding transfomer, given per unit
    S_rate = 50e3;
    Vp_rate = 2400;
    Vs_rate = 240;
    a = Vp_rate / Vs_rate;
    z_eqv_pu = 0.0217 + 1j * 0.0434; % per unit
    
    % it supplies a 24 kW load @ pf = 1, Vs is at rated value
    % input is 24.5 kW
    z_pu = Vp_rate^2 / S_rate;
    z_eqv_p = z_pu * z_eqv_pu;
    f_mdri("z_eqv_p", z_eqv_p, 1);
    
    % find Vp
    Ps = 24e3;
    pf = 1;
    tht = acos(pf);
    I_s = Ps / (Vs_rate * cos(tht));
    I_p = I_s / a;
    V_p = a*Vs_rate + I_p * z_eqv_p;
    f_mdri("V_p", V_p, 1);
    
    % efficiency...given by values
    P_out = I_s * Vs_rate; % given as 24 kW
    P_in = 24.5e3;
    eff_eta = P_out / P_in; % will reduce, need a core loss
    fprintf("efficency:  %0.1f %%\n", 100*eff_eta);
    
    % find the core copper loss, input to output
    P_loss = P_in - P_out;
    P_copper = P_loss - I_p^2 * real(z_eqv_p);
    P_core = P_loss - P_copper;
    fprintf("loss, copper:  %0.2f W\n", P_copper);
    fprintf("loss, core:  %0.2f W\n", P_core);
    
    % estimate core loss resistance, primary
    R_core = Vp_rate^2 / P_core;
    fprintf("the core resistance is:  %0.1f k Ohm\n", R_core/1e3);
    
    % estimate X_m if i_ex given..
    i_ex = 0.3;
    S_core = Vp_rate * i_ex;
    Q_core = sqrt(S_core^2 - P_core^2);
    X_m = Vp_rate^2 / Q_core;
    fprintf("Xm = %0.1f k Ohm\n", X_m / 1e3);
    
    % max eff should be used as given, because P_core = P_copper there...
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END> 
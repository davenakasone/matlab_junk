%{
    1  :  quiz 1, basic power
    2  :  quiz 2, mag circuit
    3  :  quiz 3, transformer
    4  :  quiz 4, transmission line
%}
format compact;
close all;
clear all;
clc;
select = 4;


%-------------------------------------------------------------------------------------
if select == 1
    % 120 V, 60 Hz circuit has 15 A breaker
    % powers 200 W bulb and 500 VA refrigerator  @ 80% pf, lag
    V_in = 120;
    f = 60;
    w = 2 * pi * f;
    
    pf_bulb = 1;
    theta_bulb = acos(pf_bulb);
    S_abs_bulb = 200;
    P_bulb = S_abs_bulb * cos(theta_bulb);
    Q_bulb = S_abs_bulb * sin(theta_bulb);
    S_bulb = P_bulb + 1j * Q_bulb;
    
    pf_frige = 0.8;
    theta_frige = -1 * acos(pf_frige); % lag
    S_abs_frige = 500;
    P_frige = S_abs_frige * cos(theta_frige);
    Q_frige = S_abs_frige * sin(theta_frige);
    S_frige = P_frige + 1j * Q_frige;
    
    P_total = P_bulb + P_frige;
    Q_total = Q_bulb + Q_frige;
    S_total = P_total + 1j * Q_total;
    S_abs_total = abs(S_total);
    
    % find total current supplied
    I_bulb = (S_bulb / V_in);
    I_frige = conj(S_frige / V_in);
    I_total = I_bulb + I_frige;
    fprintf("\nI_abs_total:  %0.3f A\n", abs(I_total));
    
    % find overall pf
    theta_total = angle(S_total);
    fprintf("overall angle:  %0.3f deg, pf:  %0.3f lag\n",...
        rad2deg(theta_total), cos(theta_total));
    
    % a hair dryer that draws 1350 W @ 90%pf lag is added
    % will the circuit trip?
    pf_hair = 0.9;
    theta_hair = -1 * acos(pf_hair);
    P_hair = 1350;
    S_abs_hair = P_hair / pf_hair;
    Q_hair = S_abs_hair * sin(theta_hair);
    S_hair = P_hair + 1j * Q_hair;
    S_total = S_total + S_hair;
    I_total = conj(S_total / V_in);
    f_mdri("I_total", I_total, 1);
    fprintf("I_abs =  %0.3f A , it will trip\n", abs(I_total));
    f_mdri("S_total", S_total, 1);
    I_result = real(S_total) / V_in;
    fprintf(" even with %0.3f cap, I = %0.3f, still trips\n",...
        imag(S_total), I_result);
    
    % new source current if bulb off
    S_total = S_total - S_bulb;
    I_total = conj(S_total / V_in);
    f_mdri("S_new", S_total, 1);
    f_mdri("I_new", I_total, 1);
    fprintf("it is under 15A, no trip\n");
end


%-------------------------------------------------------------------------------------
if select == 2
    % mag circuit, sheet steel
    V_peak = 170;
    V_rms = 120;
    w = 2 * pi * 60;
    
    depth_m = 15/100;
    l_m = (15/100 + 50/100) * 4;
    a_m = depth_m^2;
    n_1 = 20; % powered
    n_2 = 40; % not energized
    flux_peak = V_peak / (n_1 * w); % find peak flux
    B = flux_peak / a_m; % using phi = B A
    fprintf("\nB should be %0.3f T\n", B);
    
    % determine peak value of current supplied by source
    H = 180; % from curve
    mu = B / H; % using B = mu H
    I = H * l_m / n_1; % using NI=Hl, amper's law
    fprintf("I_peak supplied = %0.3f A\n", I);
    
    % determine peak value in non-energize coil
    E_peak = n_2 * w * flux_peak; % must share same flux
    fprintf("peak V in non-enegerized:  %0.3f\n", E_peak);
    
    % voltage supply increases 20%, current increases by what?
    % flux density is proportional to voltage
    % flux intensity is proportional to current
    flux_peak = 1.2 * V_peak / (n_1 * w);
    B = flux_peak / a_m;
    H = 300; % update from curve
    II = H * l_m / n_1;
    fprintf("the current went up: %0.3f %%   not 100?\n", 100*(II-I)/I);
    
    % flux density if AC source replaced by 170V dc,
    % internal resistance of 0.4 ohms
    % inductance is short in DC steady state
    I = 170 / 0.4; % ohm's law
    H = n_1 * I / l_m; % go to the curve
    B = 1.56;
    fprintf("flux density updated, B = %0.3f T\n", B);
end


%-------------------------------------------------------------------------------------
if select == 3
    % 2-winding transformer
    % S = 5 kVA rating
    % f = 60 Hz,  480/120 V , P_core = 35 W, i_ex = 0.22A, given z_eqv_p
    f = 60;
    w = 2 * pi * f;
    S_rated = 5e3;
    V_p = 480;
    V_s = 120;
    a = V_p / V_s;
    i_ex = 0.22;
    P_core = 35;
    z_eqv_p = 1.6 + 1j * 3.5;
    
    % use the primary side model, find R_core
    % use open circuit test
    R_core = V_p^2 / P_core; % because P = V^2 / R
    S_abs_core = V_p * i_ex; % because S = VI
    Q_core = sqrt(S_abs_core^2 - P_core^2);
    X_mag = V_p^2 / Q_core;
    fprintf("\nR_core:  %0.3f ohms,  X_m:  %0.3f\n", R_core, X_mag);
    
    % perform short circuit test, what what abs(V_p) is needed,
    % so I_p gets to rated value? ignore shunt
    I_rated = S_rated / V_p;
    V_apply = I_rated * abs(z_eqv_p);
    fprintf("apply voltage:  %0.3f V\n", V_apply);
    
    % ignore shunt, now someone shorts while connected
    % find the current
    I_short = V_p / abs(z_eqv_p);
    fprintf("I_short = %0.3f, about %d x bigger\n", I_short, I_short/I_rated);
    
    % load connected on secondary side
    % it appears as a^2 * Z_L
    % determine magnitude of V_s
    z_load = 4 + 1j*3;
    temp = abs(z_load) / abs(a^2 * z_load + z_eqv_p);
    a_Vs = V_p * a^2 * temp; % voltage division
    new_Vs = a_Vs / a;
    fprintf("Vs = %0.3f\n", new_Vs);
    
    % compute efficieny given above conditions
    I_s = new_Vs / z_load;
    P_out = real(new_Vs * conj(I_s));
    P_copper = abs(I_s/a)^2 * real(z_eqv_p);
    eta_eff = P_out / (P_out + P_core + P_copper);
    fprintf("eff = %0.3f %%\n", 100*eta_eff);
end


%-------------------------------------------------------------------------------------
if select == 4
    %{
        usually parameters to mid-length 2 port are given
            if not, call f_line_med_ABCD()
        use these identities:
            V_s = A * V_r + B * I_r
            I_s = C * V_r + D * I_r
        ...if it is short, V_r = 0 -->  V_s = B * I_r  ,  I_s = D * I_r
        ...if it is open , I_r = 0 -->  V_s = A * V_r  ,  I_s = C * V_r
        fixed sender, make Thevinin:
            V_th = V_s / A # open
            Z_th = f_para(z_line, (2/y_shunt)) # short
    %}
    % a balanced 3-phase transmission line, sender is variable, receive is fixed
    f = 60; % Hz
    w = 2 * pi * f; % rad/sec
    z_line = 1j*175; % ohms
    y_shunt = 1j*(5e-4); % S
    A = 0.956 * exp(1j*deg2rad(0));
    B = 175 * exp(1j*deg2rad(90));
    C = 4.89e-4 * exp(1j*deg2rad(90));
    D = A;
    V_s = 141.45e3 * exp(1j*deg2rad(0)); % V
    V_ll = 245e3; % the line to line voltage of sender is fixed
    z_th = f_para(z_line, 2/y_shunt); % also given
    
    % q1, use the 2 port equation
    V_r = V_s / A; % this is V_th
    fprintf("q1, Vr under no load:  %0.3f kV\n", V_r/1e3);
    
    % q2, find Is, line is shorted at the receiving end
    % I_s = D * I_r
    % V_s = B * I_r
    I_s = (V_s / B) * D;
    fprintf("\nq2, find Is if receving end is short, I_s =  %0.3f A\n", abs(I_s));
    
    % q3, find purley reactive load to make Vr = Vs = 141.45 kV
    % Vr = Vth * X/(X + z_th)
    % Vr/Vth = A     ...vth = vs/A = vr  if ==, vr/vth = A
    syms X;
    eqn_x = X / (X + imag(z_th));
    x_load = solve(eqn_x==A, X);
    fprintf("\nq3, use X =  %0.3f ohms\n", x_load);
    
    % q4, find purley resitive element, vs=vr=141.45 kV
    syms R;
    eqn_r = R^2 / (R^2 + abs(z_th)^2);
    r_load = solve(eqn_r == A^2, R);
    fprintf("\nq4, use R =  %0.3f ohms\n", r_load(2));
  
    % q5, the 2-phase max power that can be delivered to a purley resistive load
    % Vr limited to 0.95 and 1.05 pu
    Vr_95 = 0.95 * V_s;
    vth = V_s/A;
    rate = Vr_95 / vth;
    eqn = R^2 / (R^2 + imag(z_th)^2);
    r_load = solve(eqn == rate^2, R);
    use_R = r_load(2);
    p_tot = 3 * Vr_95^2 / use_R;
    fprintf("\nq5,  total power:  %0.3f  MW\n", p_tot/1e6);
    
    
    
 
    
  
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END> 
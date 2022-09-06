%{
    1  :  s22, p1, mag circuit
    2  :  s22, p2, transmission line   
            ...q7 was a mistake, Vs not Vr
            ...q12 looks fucked up, check work
    3  :  s22, p3, transformer
    
%}
format compact;
close all;
clear;
clc;
select = 3;


%-------------------------------------------------------------------------------------
if select == 1
    % core with an air-gap, assume it is linear, no loss, no fringe
    f = 60; % Hz
    w = 2 * pi * f; % rad/sec
    V_in = 240; % V RMS
    N = 200; % turns
    l_core = 1.195; % m
    l_air = 0.5/100; % m
    a_m = 30e-4; % m^2
    mu_core = 1.2e-3;
    mu_air = f_perm_mu(1);
    
    reluct = (l_core / (mu_core * a_m)) + (l_air / (mu_air * a_m));
    phi_rms = V_in / (w * N);
    
    % q1, find B_rms in the core:
    B_rms = phi_rms / a_m; % because phi = BA
    fprintf("\nq1,  B_rms =  %0.2f T\n", B_rms);
    
    % q2, find I_rms, the source current:
    I_rms = phi_rms * reluct / N; % because NI = phi * reluctance
    fprintf("\nq2,  I_rms = %0.2f A\n", I_rms);
    
    % q3, calculate the inductance reactance of the coil:
    X_coil = V_in / I_rms;
    fprintf("\nq3,  X_coil = %0.2f ohms\n", X_coil);
    
    % q4, what should N be to draw 2.4 kVAR of reactive power from the soucre:
    Q = 2.4e3;
    I = Q / V_in;
    N_new = phi_rms * reluct / I;
    fprintf("\nq4,  N = %0.2f turns\n", N_new);
    
    % q5, repeat q2, core is solid, N=200...find current:
    l_new = 1.2;
    reluct_new = l_new / (mu_core * a_m);
    I_new = phi_rms * reluct_new / N;
    fprintf("\nq5,  I = %0.2f A\n", I_new);
    
    % q6, repeat q4, core is solid...find turns:
    N_new = phi_rms * reluct_new / I;
    fprintf("\nq6,  N = %0.2f turns\n", N_new);
end


%-------------------------------------------------------------------------------------
if select == 2
    % a balanced 3-phase transmission line, sender is variable, receive is fixed
    f = 60; % Hz
    w = 2 * pi * f; % rad/sec
    z_line = 5 + 1j*5; % ohms
    y_shunt = 1j*0.0008; % S
    A = 0.98 * exp(1j*deg2rad(0));
    B = 50.25 * exp(1j*deg2rad(84.3));
    C = 0.00079 * exp(1j*deg2rad(90));
    D = A;
    V_r = 130e3 * exp(1j*deg2rad(0)); % V
    V_ll = 225e3; % the line to line voltage...
    
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
    
    % q7, find V_s under no load
    V_s = A * V_r;
    fprintf("\nq7, V_s =  %0.2f kV\n", V_s/1000);
    
    % q8, find I_s...sending current under no load / "charging"
    I_s = C * V_r;
    fprintf("\nq8, I_s =  %0.2f A\n", imag(I_s)); % pure imag...
    
    % a 3-phase load draws 240 MW and 90 MVAR, on the receiving end
    P_tot = 240e6;
    Q_tot = 90e6;
    S_tot = P_tot + 1j*Q_tot;
    fprintf("\n");
    f_mdri("S_tot", S_tot, 1/10^6);
    
    % q9, magnitude and phase andgle of current at receiving end
    I_r = conj(S_tot / (sqrt(3) * V_ll)); % must conjugate...
    fprintf("\nq9, ");
    f_mdri("I_r", I_r, 1);
    
    % q10, find the per-phase voltage at the sending end
    V_s = A * V_r + B * I_r;
    fprintf("\nq10, ");
    f_mdri("V_s", V_s, 1/1000);
    
    % q11, find the mag and phase of sending-end current
    I_s = C * V_r + D * I_r;
    fprintf("\nq11, ");
    f_mdri("I_s", I_s, 1);
    
    % q12, find the total power supplied to line at sending end
    S_in = 3 * V_s * conj(I_s);
    P_in = 3 * abs(V_s) * abs(I_s) * cos(angle(V_s) - angle(I_s));
    Q_in = 3 * abs(V_s) * abs(I_s) * sin(angle(V_s) - angle(I_s));
    check_P = real(S_in) - P_in
    check_Q = imag(S_in) - Q_in
    fprintf("\nq12, P =  %0.1f W  ,  Q =  %0.1f MVAR\n",...
        P_in/10^6, Q_in/10^6);
end


%-------------------------------------------------------------------------------------
if select == 3
    % single phase, 2 winding, shunt unknown
    S_rate = 50e3; % VA
    Vp_rate = 2400; % primary voltage
    Vs_rate = 240;
    a = Vp_rate / Vs_rate; % turns ratio
    z_series_pu = 0.0217 + 1j*0.434; % per unit 
    P_in = 24.4e3; % W
    P_out = 24e3; % W, while suppling load with pf = 1
    
    % q13, find z_series on the primary side
    z_base = Vp_rate^2 / S_rate;
    z_series = z_base * z_series_pu;
    fprintf("q13, ");
    f_mdri("z_series", z_series, 1);
    
    % q14, find the primary side voltage, KVL
    I_p = P_out / Vp_rate;
    V_p = Vp_rate + z_series * I_p;
    fprintf("\nq14, ");
    f_mdri("VP", V_p, 1);
    
    % q15, determine core losses
    P_loss = P_in - P_out;
    P_copper = a^2 * real(z_series);
    P_core = P_loss - P_copper; 
    fprintf("\nq15, core losses:  %0.3f  W\n", P_core);
    
    % q16, estimate core loss resistance, primary side
    R_core = Vp_rate^2 / P_core;
    fprintf("\nq16 , core resistance:  %0.3f k ohms\n", R_core/1e3);
    i_core = P_core / Vp_rate;
    
    % q17, estimate transform i_e, no load + mag reactance X_m = 20 k
    X_m = 20e3;
    i_m = Vp_rate / X_m;
    i_ex = sqrt(i_m^2 + i_core^2);
    fprintf("\nq17 , i_m =  %0.3f  mA\n", i_ex*1e3);
    
    % q18, find VR and eff
    VR = 100 * (real(V_p) - Vp_rate) / Vp_rate;
    eff_eta = 100 * (P_out/P_in);
    fprintf("\nq18 , VR:  %0.3f %%  ,  eff:  %0.3f %%\n", VR, eff_eta);
    
end


%-------------------------------------------------------------------------------------
if select == 99
    
end
%%%%%%%%~~~~~~~~~END> 
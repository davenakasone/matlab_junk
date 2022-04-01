%{
    1  :  slide9, ex1, ideal transformer
    2  :  slide25, ex2, real transformer, given oc and sc test results
    3  :  slide27, ex3, make base units 
    4  :  slide32, ex4, applications ???
    5  :  p3.1, a transformer
    6  :  p3.2, model of high side given
    7  :  p3.3, open circuit, closed circuit info is given
    8  :  p3.4, practical transformer
    9  :  p3.5, different power in countries
    10 :  p3.6, transformer, basic
    11 :  p3.7, transformer, per-unit
    12 :  p3.8, transformer, per-unit
    13 :  hw3, the transformer design
    14 :  quiz 3   
%}
close all;
clear all;
clc;
select = 13;


%-------------------------------------------------------------------------------------
if select == 1
    z_line = 0.18 + 1j*0.24;
    z_load = 4 + 1j*3;
    v_p = 480*exp(1j*deg2rad(0));
    i_g = v_p / ( z_line + z_load);
    f_mdri("I_G", i_g, 1);
    v_load = i_g * z_load;
    f_mdri("V_load", v_load, 1);
    p_loss = abs(i_g)^2 * real(z_line); % peak current ^ 2  * resistance
    f_mdri("P_loss", p_loss, 1);
    
    % a 1:10 step-up on generator, 10:1 step down on load
    a_2 = 10/1;
    a_1 = 1/10;
    z_load_p = z_load * a_2^2; % eliminate T2, reflect impedance
    z_sum = z_line + z_load_p; % add series impedance
    z_sum_p = z_sum * a_1^2; % eliminate T1, reflect impedance to source
    i_p = v_p / z_sum_p;
    f_mdri("I_G", i_p, 1);
    
    % now that you know source current and the 2 turns ratios, get line and load:
    i_line = a_1 * i_p;
    f_mdri("I_line", i_line, 1);
    i_load = a_2 * i_line;
    f_mdri("I_load", i_load, 1);
    % the load voltage is achieved:
    v_load = i_load * z_load;
    f_mdri("V_load", v_load, 1);
    p_loss = abs(i_line)^2 * real(z_line);
    f_mdri("P_loss", p_loss, 1);
end


%-------------------------------------------------------------------------------------
if select == 2
    S = 20e3; % VA
    f = 60; % Hz
    w = 2 * pi * f; % rad/sec
    V_oc = 8000;
    I_oc = 0.214;
    P_oc = 400;
    V_sc = 489;
    I_sc = 2.5;
    P_sc = 240;
    a = 8000/240; % turns ratio
end


%-------------------------------------------------------------------------------------
if select == 3
    a = 8000/240;
    f = 60;
    w = 2 * pi * f;
    S = 20e3;
    R_c = 159e3;
    X_m = 38.4e3;
    R_eq = 38.3;
    X_eq = 192;
    
    V_base = 8000;
    S_base = 20e3;
    Z_base = V_base^2 / S_base;
    Z_SE_pu = (R_eq + 1j*X_eq) / Z_base; % per unit
    R_c_pu = R_c / Z_base; % per unit
    X_m_pu = X_m / Z_base; % per unit
end


%-------------------------------------------------------------------------------------
if select == 4
    S = 15e3;
    a = 2300/230;
    V_oc = 2300;
    I_oc = 0.21;
    P_oc = 50;
    V_sc = 47;
    I_sc = 6;
    P_sc = 160;
    
    % finding equiv for high side
    % finding equiv for low side
    % calculate pfs
    % efficiency
end


%-------------------------------------------------------------------------------------
if select == 5
    a = 50/200;
    f = 60;
    w = 2 * pi * f;
    v_s = (282.8/sqrt(2)) * exp(1j*deg2rad(0));
    f_mdri("V_s rms", v_s, 1);
    v_sp = a * v_s;
    i_s = (7.07/sqrt(2)) * exp(1j*deg2rad(-36.87));
    f_mdri("I_s rms", i_s, 1);
    i_sp = i_s / a;
    R_eq = 0.05;
    R_c = 75;
    X_eq = 0.225;
    X_m = 20;
    
    % apply KVl for Vp
    v_p = v_sp + i_sp * (R_eq + 1j*X_eq);
    f_mdri("V_p", v_p, 1);
    
    % apply KCL for Iex
    i_ex = (v_p / R_c) + (v_p / (1j*X_m));
    f_mdri("I_ex", i_ex, 1);
    
    % apply KCL for Ip
    i_p = i_ex + i_sp;
    f_mdri("I_p", i_p, 1);
    
    % find VR at this load
    VR = (abs(v_p) - abs(a * v_s)) / abs(a * v_s);
    fprintf("\nVR  %0.2f  percent\n", 100*VR);
    
    %input power
    S_in = v_p * conj(i_p); % another option...take real part
    P_in = abs(v_p) * abs(i_p) * cos(angle(v_p) - angle(i_p));
    f_mdri("P_in", P_in, 1);

    %output power
    S_out = v_s * conj(i_s); % total solution in S
    P_out = abs(v_s) * abs(i_s) * cos(angle(v_s) - angle(i_s));
    f_mdri("P_out", P_out, 1);
    
    % efficiency
    eff_eta = P_out / P_in
end


%-------------------------------------------------------------------------------------
if select == 6
    S_base = 20e3;
    V_base = 8000;
    a = 8000/277;
    R_p = 32;
    X_p = 45;
    R_c = 250e3;
    R_s = 0.05;
    X_s = 0.06;
    X_m = 30e3;
    
    % make equiv circuit on high side
    R_s_p = a^2 * R_s
    X_s_p = a^2 * X_s
    
    % find the per-unit equiv...use Vp and S given
    I_base = S_base / V_base
    Z_base = V_base / I_base
    R_s_pu = R_s_p / Z_base % uses equiv value reflected
    X_s_pu = X_s_p / Z_base % uses equiv value reflected
    X_m_pu = X_m / Z_base
    R_p_pu = R_p / Z_base
    X_p_pu = X_p / Z_base
    R_c_pu = R_c / Z_base
    
    % transformer supplies load 277 V @ 0.8 pf, lag...find input voltage
    V_s_mag = 277; % it operates at rating
    pf = 0.8; % remember it is lagging, so it is negative    th_v - th_i < 0
    theta = -1 * acos(pf);
    V_s = V_s_mag * exp(1j*0); % reference at 0
    V_s_p = a * V_s;
    f_mdri("V_s_p", V_s_p, 1);
    I_s_mag = S_base / V_s_mag;
    I_s = I_s_mag * exp(1j*theta);
    f_mdri("I_s", I_s, 1);
    I_s_p = I_s / a;
    f_mdri("I_s_p", I_s_p, 1);
    R_eq = R_p + R_s_p;
    X_eq = 1j*(X_p + X_s_p);
    V_p = V_s_p + I_s_p * (R_eq + X_eq); % simplified KVL, ignores shunt
    f_mdri("V_p", V_p, 1);
    VR = 100 * (abs(V_p) - V_base) / V_base % percent
    
    % under above conditions, find losses   ...lots of assumptions
    P_out = S_base * cos(theta)
    P_out = real(V_s_p * conj(I_s_p)) % either
    P_copper = abs(I_s_p)^2 * R_eq % assumes same current hits both resistors
    P_core = abs(V_p)^2 / R_c % core uses derived primary voltage
    eff_etta = P_out / (P_out + P_copper + P_core) 
end


%-------------------------------------------------------------------------------------
if select == 7
    S = 2e3; % all data is given on primary side
    a = 230/115;
    V_oc = 230;
    I_oc = 0.45;
    P_oc = 30;
    V_sc = 13.2;
    I_sc = 6;
    P_sc = 20.1;
    
    % using open circut, shunt is analyzed ...solve through Y
    Y_ex_mag = I_oc / V_oc;
    theta = -1 * acos(P_oc / (V_oc * I_oc)); % for invert
    Y_ex = Y_ex_mag * exp(1j*theta);
    f_mdri("Y_ex", Y_ex, 1);
    R_c = 1 / real(Y_ex);
    X_m = 1 / (-1 * imag(Y_ex)); % for invert
    
    % using short circuit, line components solved through Z
    Z_eq_mag = V_sc / I_sc;
    theta = acos(P_sc / (V_sc * I_sc));
    Z_eq = Z_eq_mag * exp(1j*theta);
    f_mdri("Z_eq", Z_eq, 1);
    R_eq = real(Z_eq);
    X_eq = imag(Z_eq);
    
    % convert circuit to secondary side
    R_eq_s = R_eq / a^2;
    X_eq_s = X_eq / a^2;
    R_c_s = R_c / a^2;
    X_m_s = X_m / a^2;
    Z_eq_s = R_eq_s + 1j*X_eq_s;
    
    S_s = S / a; % it is on secondary side
    V_p = 230; % given
    V_s = V_p / a; % given
    V_p_s = V_p / a; % on secondary side
    I_s_mag = S_s / V_p_s; % solving S = VI
    %VR @ 1 unity
    pf = 1;
    theta = acos(pf); % unity
    I_s = I_s_mag * exp(1j*theta);
    V_ps = V_s + Z_eq_s * I_s;
    f_mdri("V_ps", V_ps, 1);
    VR = ((abs(V_ps) - V_s) / V_s) * 100 % percent
    %VR @ 0.8 lead
    pf = 0.8;
    theta = acos(pf); % positive for lead
    I_s = I_s_mag * exp(1j*theta);
    V_ps = V_s + Z_eq_s * I_s;
    f_mdri("V_ps", V_ps, 1);
    VR = ((abs(V_ps) - V_s) / V_s) * 100 % percent
    %VR @ 0.8 lag, by KVL
    pf = 0.8;
    theta = -1 * acos(pf); % negative for lag
    I_s = I_s_mag * exp(1j*theta);
    V_ps = V_s + Z_eq_s * I_s;
    f_mdri("V_ps", V_ps, 1);
    VR = ((abs(V_ps) - V_s) / V_s) * 100 % percent
    
    %efficiency at 0.8 pf lag
    pf = 0.8;
    theta = -1 * acos(pf);
    P_out = V_s * I_s_mag * cos(theta)
    P_copper = I_s_mag^2 * R_eq_s
    P_core = abs(V_ps)^2 / R_c_s % see above
    eff_eta = P_out / (P_out + P_copper + P_core)
end


%-------------------------------------------------------------------------------------
if select == 8
    S = 100e3;
    V_p = 14e3;
    V_s = 2.4e3;
    a = V_p / V_s;
    R_p = 38.2;
    X_p = 140;
    R_s = 0.12;
    X_s = 0.5;
    P_load = 90e3;
    pf_load = 0.9; % lag == negative angle
    theta_load = -1 * acos(pf_load);
    V_load_mag = 2300;
    
    % consolidate given info
    Z_p = R_p + 1j*X_p;
    Z_p_s = Z_p / a^2;
    Z_s = R_s + 1j*X_s;
    Z_eq = Z_p_s + Z_s;

    I_s_mag = P_load / (V_load_mag * pf_load);
    I_s = I_s_mag * exp(1j*theta_load);
    f_mdri("I_s", I_s, 1);
    
    % find voltage on secondary side, KVL low side circut
    V_ps = V_load_mag  + I_s * Z_eq;
    f_mdri("V_ps", V_ps, 1); % primary voltage on secondary side
    V_p = a * V_ps;
    f_mdri("V_p", V_p, 1); % voltage at the source
    
    % to find the VR, need Vp on primary under full load, KVL
    V_ps = V_load_mag + I_s * Z_s;
    f_mdri("V_ps", V_ps, 1);
    VR = (abs(V_ps) - V_load_mag) / V_load_mag
    
    % efficiency by P_out/P_in
    P_out = P_load % given
    V_ps = V_load_mag  + I_s * Z_eq;
    P_in = real(V_ps * conj(I_s))
    eff_eta = P_out / P_in
end


%-------------------------------------------------------------------------------------
if select == 9
    na_Vrms = 120; % wall voltage in north america
    na_f = 60;
    na_w = 2 * pi * na_f;
    eu_Vrms = 230; % wall voltage in EU is 220-240
    eu_f = 50;
    eu_w = 2 * pi * eu_f;
    S = 1000;
    N_120 = 500; % turns on 120V side
    N_240 = 1000; % turns on 240V side
    
    % when supplying through 120 V, and no 240 V load
    % flux in core is:  phi(t) = -V_m / ( w * N_120)   cos(wt)
    % use peaks, not RMS
end


%-------------------------------------------------------------------------------------
if select == 10
    S = 15e3;
    Vp = 8000;
    Vs = 230;
    a = Vp / Vs;
    Z_p = 80 + 1j*300;
    R_c = 350e3;
    X_m = 70e3;
    Z_ex = R_c + 1j* X_m;
    
    % calculate for part A, reflect Z_load to primary side
    V_p = 7967; % given
    Z_load = 3.2 + 1j*1.5;
    Z_load_p = a^2 * Z_load;
    I_s_p = V_p / (Z_load_p + Z_p);
    f_mdri("I_s_p", I_s_p, 1);
    V_s_p = I_s_p * Z_load_p;
    f_mdri("V_s_p", V_s_p, 1);
    V_s = V_s_p / a;
    f_mdri("V_s", V_s, 1);
    VR = 100*(V_p - abs(V_s_p)) / abs(V_s_p)
    
    % calculate for part B, reflect Z_load to primary side
    V_p = 7967; % given
    Z_load = 0 - 1j*3.5;
    Z_load_p = a^2 * Z_load;
    I_s_p = V_p / (Z_load_p + Z_p);
    f_mdri("I_s_p", I_s_p, 1);
    V_s_p = I_s_p * Z_load_p;
    f_mdri("V_s_p", V_s_p, 1);
    V_s = V_s_p / a;
    f_mdri("V_s", V_s, 1);
    VR = 100*(V_p - abs(V_s_p)) / abs(V_s_p)
end


%-------------------------------------------------------------------------------------
if select == 11
    R_pu = 0.01;
    X_pu = 0.05;
    S_base = 5000e3; % given
    a = 230/13.8e3;
    % open circuit on low-side
    V_oc = 13.8e3;
    I_oc = 15.1;
    P_oc = 44.9e3;
    
    % use OC on low side to find components
    Y_ex_mag = I_oc / V_oc;
    theta = -1 * acos(P_oc / (V_oc * I_oc)); % for inversion
    Y_ex = Y_ex_mag * exp(1j*theta);
    f_mdri("Y_ex", Y_ex, 1);
    R_c = 1 / real(Y_ex)
    X_m = -1 / imag(Y_ex)
    V_base = V_oc; % infered
    Z_base = V_base^2 / S_base
    R_eq = R_pu * Z_base
    X_eq = X_pu * Z_base
    
    % load is 4000kW, 0.8pf lag, Vs= 13.8 kV,
    P_load = 4000e3;
    Vs = 13.8e3;
    pf = 0.8;
    theta = -1 * acos(pf);
    Is_mag = P_load / (Vs * pf);
    Is = Is_mag * exp(1j*theta);
    f_mdri("Is", Is, 1);
    Vp_s = Vs + Is * (R_eq + 1j*X_eq);
    f_mdri("Vp_s", Vp_s, 1);
    VR = (abs(Vp_s) - Vs) / Vs
    P_copper = abs(Is)^2 * R_eq
    P_core = abs(Vp_s)^2 / R_c
    eff_eta = P_load / (P_load + P_copper + P_core)
end


%-------------------------------------------------------------------------------------
if select == 12
    S_base = 150e6;
    V_lo = 15e3;
    V_hi = 200e3;
    a = V_lo / V_hi;
    R_pu = 0.012;
    X_pu = 0.05;
    X_mag_pu = 100;
    
    % the low-voltage side is the "primary"
    Z_base = V_lo^2 / S_base;
    R_eq_pu = R_pu * Z_base;
    X_eq_pu = X_pu * Z_base;
    X_m_pu = X_mag_pu * Z_base;
    Z_eq_pu = R_eq_pu + 1j*X_eq_pu;
    
    % find VR for full load current @ 0.8 lag
    % that means secondary/hi gets 150 MVA
    pf = 0.8;
    theta = -1 * acos(pf);
    I_s_p_mag = S_base / (V_lo * pf);
    I_s_p = I_s_p_mag * exp(1j*theta);
    f_mdri("I_s_p", I_s_p, 1);
    V_p = V_lo + I_s_p * Z_eq_pu;
    f_mdri("V_p", V_p, 1);
    VR = (abs(V_p) - V_lo) / V_lo
end

%-------------------------------------------------------------------------------------
if select == 13
    f = 60;             % Hz
    w = 2 * pi * f;     % rad/sec
    S = 5e3;            % VA
    V_hi = 480;         % V on high side
    V_lo = 120;         % V on low side
    a = V_hi / V_lo;    % appears to be a single phase, step-down
                        % sheet steel core, see B-H curve,
                            % peak value of operating flux density B == 1 Tesla
    P_core = 35;        % W, assumed core loss
                        % i_ex < 4%     i_ex = i_mag + i_he
                        % efficiency at full-load, unity pf  > 95.5%
   
    % determine: core size and turns; limit weight
    % use a uniform and square core + equal depth
    syms x;                  % depth, face width
    syms y;                  % inside square dimensions
    cross_area = x^2;        % cross-sectional area
    l_flux = 4 * (x + y);    % l = average length flux travles
    B = 1;                   % T, from B-H curve for sheet steel
    H = 180;                 % At/m, from B-H curve for sheet steel
    
    % preliminary analysis 1
    a = V_hi / V_lo;    % == N_p / N_s, must have a = 4
    syms N_s;           % variable for secondary turns
    N_p = a * N_s;      % must be true to maintian turns ratio
    N_t = N_p + N_s;    % total turns, both coils
    
    % preliminary analysis 2, flux to peak voltage
    % phi = B / A = B / x^2 = V / (w N_p) = V / (w * a * N_s)
    xx = sqrt(0.45/N_s); % x = sqrt(0.45/N_s), eqn1, this looks wrong
    
    % preliminary analysis 3, wires and currents
    I_p_abs = S / V_hi; % 10.4167 A, rated voltage primary side
    I_s_abs = 4 * I_p_abs; % secondary should handle 4x more: 41.6667 A
    % using AWG #8, for simplicity
    R_awg8 = 2.06e-3; % omhs / m
    A_awg8 = 8.36e-6; % m^2
    
    % preliminary analysis 4, RMS at excitation
    i_core = sqrt(2) * (P_core / V_hi); % 0.1031 A
    i_mag = H * (x + y) / N_s; % I_m = H l / Np
    % RMS value of i_ex as % of rated current
    i_ex = 100 * sqrt( i_core^2 + i_mag^2 ) / (sqrt(2)*I_p_abs); % eqn2
    
    % preliminary analysis 5, space
    % assuming that 25% of coild space is insulation
    % total cross section area of both windings = 0.5625 y^2
    yy = sqrt(N_s / (13.457 * 10^3));    % eqn_3,  y = space constraint
    
    % preliminary analysis 6, winding resistance
    R_s = N_s * (4 * x + 0.8 * y) * R_awg8;
    R_p = 4 * N_s * (4 * x + 4.8 * y) * R_awg8;
    R_eq_p = R_p + a^2 * R_s;    % reflect to primary
    
    % preliminary analysis 7, efficiency
    % pf = 1, i_load = I_p_abs = 10.417 A
    P_out = S; % 5 kW, all real power
    P_in = P_out + P_core + R_eq_p * I_p_abs^2;
    eff_eta = 100 * (P_out / P_in); % efficency
    
    % preliminary analysis 8, find weight in kg
    t_kg = (31400 * x^2 + 1.495 * N_s) * (x + y); % tranformer weight in kg
    
    % find: x, y, Ns, Np = a Ns
    % satisfy: i_ex < 4%, eff >= 95.5% @ pf = 1
    % very Ns[1:200]
    takeNs = 200;
    Ns_vect = 1:1:takeNs;
    %
    fprintf("\n Ns    |  Np    |  x     |  y      |  i_ex/i_p  |  R_eq_p  |  eff  | weight\n");
    fprintf("turns  |turns   |meters  |meters   |%% RMS       |ohms      |%%      |kg\n");
    fprintf("--------------------------------------------------------------------------------\n");
    %for ii = 1:1:takeNs
    for ii = 85:1:105
        temp_Ns = Ns_vect(ii);
        temp_Np = subs(N_p, N_s, temp_Ns);
        temp_xx = subs(xx, N_s, temp_Ns);
        temp_yy = subs(yy, N_s, temp_Ns);
        temp_i_ex = subs(i_ex, [x, y, N_s], [temp_xx, temp_yy, temp_Ns]);
        temp_R_eq_p = subs(R_eq_p, [x, y, N_s], [temp_xx, temp_yy, temp_Ns]);
        temp_eff_eta = subs(eff_eta, [x, y, N_s], [temp_xx, temp_yy, temp_Ns]);
        temp_t_kg = subs(t_kg, [x, y, N_s], [temp_xx, temp_yy, temp_Ns]);
        fprintf("%4d   | %4d   | %6.3f | %6.3f  | %9.3f  | %6.3f   | %0.2f | %0.3f\n",...
            temp_Ns, temp_Np, temp_xx, temp_yy, temp_i_ex, temp_R_eq_p, temp_eff_eta, temp_t_kg);
    end
    %
    %
    select_Ns = 93; % select optimal, based on table
    new_Np = double(subs(N_p, N_s, select_Ns));
    new_x = double(subs(xx, N_s, select_Ns));
    new_y = double(subs(yy, N_s, select_Ns));
    new_i_mag = double(subs(i_mag, [x, y, N_s], [new_x, new_y, select_Ns]));
    new_iex_rms = double(subs(i_ex, [x, y, N_s], [new_x, new_y, select_Ns])); % RMS % of max
    new_R_eq_p = double(subs(R_eq_p, [x, y, N_s], [new_x, new_y, select_Ns]));
    new_eff_eta = double(subs(eff_eta, [x, y, N_s], [new_x, new_y, select_Ns]));
    fprintf("\nmost optimal:\n");
    fprintf("Np =  %d  ,  Ns =   %d\n", new_Np, select_Ns);
    fprintf("x =  %0.3f cm  ,  y =  %0.3f cm\n", new_x * 100, new_y * 100);
    fprintf("i_ex as %% RMS:  %0.3f %%\n", new_iex_rms);
    fprintf("R_eq_p:  %0.3f  ohms\n", new_R_eq_p);
    fprintf("nominal efficiency:  %0.3f  %%\n", new_eff_eta);
    
    syms Ip;
    new_pf = 1;
    new_Pout = V_hi * Ip * new_pf;
    eqn_eta = (100 * new_Pout) / (P_core + Ip^2 * new_R_eq_p + new_Pout);
    d_eqn_eta = diff(eqn_eta, Ip);
    Ip_max = double(solve(d_eqn_eta==0, Ip>0)); % derivative == 0 --> max
    max_eta = double(subs(eqn_eta, Ip, Ip_max));
    fprintf("\nmaximum efficiency =  %0.3f %%   @  Is =  %0.3f A\n", max_eta, Ip_max);
    full_load = V_hi * I_p_abs;
    maxIs_load = V_hi * Ip_max;
    relative = (100 * maxIs_load) / full_load;
    fprintf("\nat max efficieny, %0.3f %% of the full load is used\n", relative);
    
    check_eta = (100 * Ip_max * V_hi) / (2 * P_core + V_hi * Ip_max);
    validate = check_eta - max_eta; % effectivley 0
    fprintf("\nthe maximum efficency is confirmed:  %0.3f %%\n", check_eta);
end


%-------------------------------------------------------------------------------------
if select == 14
    S = 5e3;
    f = 60;
    w = 2 * pi * f;
    V_hi = 480;
    V_lo = 120;
    a = V_hi / V_lo;
    P_core = 35;
    i_ex = 0.22;
    R_eq_p = 1.6 + 1j*3.5;
    
    R_copper = V_hi^2 / P_core
    i_he = sqrt(P_core / R_copper)
    i_m = i_ex - i_he
    
    I_p = S / V_hi
    I_sa = I_p - i_ex
   
end


%-------------------------------------------------------------------------------------
if select == 15
 
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END>  
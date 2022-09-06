%{
    1  : AC circuit using phasors
    2  : active, reactive, apparent
    3  : circuit with cable impedance
    4  : given load, use pf
    5  : power factor correction    make the triangle when in doubt
    6  : power on breakers
    7  : power consumption properties, like KCL, KVL
    8  : p8.47, ee221  solve for capacitor  ???
    9  : p8.45, symbols, unity (0 imag part) ..matlab won't help
    10 : p8.43 multi load   ???
%}
clc;
select = 10;

%-------------------------------------------------------------------------------------
if select == 1
    V_max = 100;
    phaze = deg2rad(30);
    frq = 60;
    w = 2 * pi * frq;
    V = V_max * exp(1j*phaze);
    
    R = 4;
    L = 7.95e-3;
    Z_c = 1j * w * L;
    Z = R + Z_c;
    
    Z_mag = abs(Z);
    Z_rad = angle(Z);
    Z_deg = rad2deg(Z_rad);
    
    I = V / Z;
    I_mag = abs(I)
    I_rad = angle(I);
    I_deg = rad2deg(I_rad)
    fprintf("\n\ti(t) = %f cos(wt %f)\n", I_mag, I_deg);
end

%-------------------------------------------------------------------------------------
if select == 2
    % S = VI*
    % S = P + jQ                            complex, make P, Q, S triangle
    % P = |S|cos(theta) = |S|pf             real: W
    % Q = |S|sin(theta) = |S|sqrt(1-pf^2)   reactive: var
    % |S| = |V||I| = sqrt(P^2 + Q^2)        apparent: VA
    % ....pf = P/|S| = cos(tht)    tht = acos(pf)   +/-    -lead, +lag ???
    
    % load draws 100kW, leading pf = 0.85, that means theta < 0
    pf = 0.85;                 % given
    P = 100e3;                 % given as real power = 100kW
    th_rad = -1 * acos(pf);    % must have theta < 0 if leading
    th_deg = rad2deg(th_rad); 
    S_mag = P/pf;
    Q = S_mag * sin(th_rad);
    S = P + 1j*Q;
    fprintf("\n\tP = %0.1f kW  ,  Q = %0.1f kVar  ,  |S| = %0.1f kVA\n",...
        real(S)/(1e3), imag(S)/(1e3), abs(S)/(1e3));
end

%-------------------------------------------------------------------------------------
if select == 3
    % circuit has unknow voltage source
    % 40kV ang(0) given across R(100)
    % finding I = 40k / 100 = 400 ang(0) A
    % KVL with inffered I
    R = 100;
    V_o = 40e3;
    I = V_o / R;
    Z = 5 + 1j*40;
    V_z = I * Z;
    V = V_z + V_o;
    S = V * conj(I);
    fprintf("\n\t S = %0.2f MVA  @  %0.1f deg ,  P = %0.1f MW  ,  Q = %0.1f Mvar\n",...
        abs(S)/(1e6), rad2deg(angle(S)), real(S)/(1e6), imag(S)/(1e6));
end

%-------------------------------------------------------------------------------------
if select == 4
    Z_load_R = 100;
    Z_load_L = 1j*100;
    Z_load = f_para(Z_load_R, Z_load_L); % this is a 45 degree "lag"...inductive, or positive theta
    fprintf("\n|Z_load| = %0.1f  @ %0.2f deg  ,  real= %0.1f  ,  imag= %0.1f\n",...
        abs(Z_load), rad2deg(angle(Z_load)), real(Z_load), imag(Z_load));
    
    theta = angle(Z_load);
    pf = cos(theta);
    V_o = 40e3; % @ 0 deg
    I = V_o / Z_load;
    fprintf("\n|I| = %0.1f  @ %0.2f deg  ,  real= %0.1f  ,  imag= %0.1f\n",...
        abs(I), rad2deg(angle(I)), real(I), imag(I));
    
    Z_line = 4 + 1j*40;
    V_line = I * Z_line;
    V = V_line + V_o;
    fprintf("\n|V| = %0.1f  @ %0.2f deg  ,  real= %0.1f  ,  imag= %0.1f\n",...
        abs(V), rad2deg(angle(V)), real(V), imag(V));
    S = V * conj(I);
    fprintf("\n\t|S| = %0.1f @ %0.1f deg  ,  P = %0.1f  ,  Q = %0.1f\n",...
        abs(S)/(1e6), rad2deg(angle(S)), real(S)/(1e6), imag(S)/(1e6));
end

%-------------------------------------------------------------------------------------
if select == 5
    S_mag = 100e3; % kVA
    pf = 0.8; % lagging ...that means positive
    theta_r = acos(pf); % stays positive...lagging
    theta_d = rad2deg(theta_r);
    P = S_mag * cos(theta_r); % solve using properties...atan2()
    Q = S_mag * sin(theta_r);
    
    pff = 0.95; % want to make it better, still lagging
    thetaa_r = acos(pff);
    thetaa_d = rad2deg(thetaa_r); 
    
    % Q_cap = (-1) * (P * tan(th_new) - Q)
    Q_cap = -(P * tan(thetaa_r) - Q);
    fprintf("\n\tuse a capacitor that makes %0.1f  kvar\n", Q_cap/1000);
end

%-------------------------------------------------------------------------------------
if select == 6
    % pool pump on 240 V, 60 Hz, 15A breaker {trips if > 15A}
    % motor draws 3 kW with 75% lag
    v_in = 240;
    i_lim = 15;
    f_hz = 60;
    w = 2 * pi * f_hz;
    pf = 0.75;
    P = 3e3;
    S_mag = P / pf ; 
    I_mag = S_mag / v_in;  % just had to use properties
    fprintf("I_mag = %0.1f  ,  I_lim = %0.1f ...trip\n", I_mag, i_lim);
    
    % what is equiv series impedance of motor?
    syms RR;
    R_r = double(solve(RR * I_mag^2 == P, RR)); % using power calc VI, V^2/R, or I^2 R
    fprintf("resistor on P is about  %0.1f ohms\n", R_r);
    Q = sqrt(S_mag^2 - P^2);
    R_i = double(solve(RR * I_mag^2 == Q, RR));
    fprintf("resistor in Q is about  %0.1f ohms\n", R_i);
    Z_motor = R_r + 1j*R_i;
    fprintf("|Z_motor| = %0.1f @ %0.1f deg  >>>  %0.1f + j %0.1f\n",...
        abs(Z_motor), rad2deg(angle(Z_motor)), real(Z_motor), imag(Z_motor));
    
    verify = v_in / Z_motor;
    v_pf = cos(angle(verify));
    fprintf("\nverifiying, the pf is the same:  %0.2f\n", v_pf);
    
    % shunt cap in vars, need source current to equal breaker 15A
    S_mag_new = 15 * 240; % |S| = |V||I|
    Q_new = sqrt(S_mag_new^2 - P^2); % implied new var
    S_new = P + 1j*Q_new; % get new power 
    tht_new_r = angle(S_new); % find new power's angle
    Q_cap = -(P * tan(tht_new_r) - Q); % use the formula
    pf_new = P/abs(S_new);
    fprintf("\nthe new pf:  %0.2f\n", pf_new);
    
    % get a bank rated at 1.5 k Var, what pf results?
    Q_bank = 1.5e3;
    Q_newb = Q - Q_bank;
    S_mag_new = sqrt(P^2 + Q_newb^2);
    pf_new = P / abs(S_mag_new);
    fprintf("\nbank makes pf=  %0.2f lag\n", pf_new);
    
    % optimal, Q = 0
    Q_new = 0;  % ...Q - Qnew = -Qnew, C = -Qc / wV^2
    C_opt = Q / (w * v_in^2);
    fprintf("\nmake cap  %0.1f uF\n", C_opt*(1e6));
    
end

%-------------------------------------------------------------------------------------
if select == 7
    % ciruit has unknonw source, 
    % given load voltage, line impedance, and load power
    Z_line = 0.1 + 1j*0.2;
    v_load = 115;
    S_mag_load = 2.3e3;
    pf_load = 0.8; % lag
    
    % real and reactive power consumed by load
    tht_r = acos(pf_load);
    P_load = S_mag_load * cos(tht_r);
    fprintf("\nload real pow, P =  %0.1f W\n", P_load);
    P_load = S_mag_load * pf_load; % same as above
    fprintf("load real pow, P =  %0.1f W\n", P_load);
    Q_load = S_mag_load * sqrt(1-pf_load^2);
    fprintf("load reac pow, W =  %0.1f W\n", Q_load);
    Q_load = S_mag_load * sin(tht_r); % same as above
    fprintf("load reac pow, W =  %0.1f W\n", Q_load);
    S_load = P_load + 1j*Q_load;
    check = abs(S_load);
    checkk = cos(angle(S_load));
    I_mag = abs(S_load) / v_load;
    Z_load = (P_load/I_mag^2) + 1j*(Q_load/I_mag^2);
    V_load = v_load * exp(1j*tht_r);
    f_mdri("V_load", V_load, 1/1000);
    
    % real and reactive power consumed by line
    % ... must have same current
    V_mag = I_mag * abs(Z_line);
    P_line = I_mag^2 * real(Z_line); % P = I^2 R
    Q_line = I_mag^2 * imag(Z_line); % Q = I^2 Z_L
    S_line = P_line + 1j*Q_line;
    fprintf("\nin line, P = %0.1f W, Q = %0.1f var\n", P_line, Q_line);
    
    % and solving through the source...just conserve power
    S_source = S_load + S_line; % technically negative since a supply of power
    V_s_mag = abs(S_source) / I_mag;
    V_s = I_mag * (Z_line + Z_load);
    f_mdri("V_s", V_s, 1); 
    tht_diff = rad2deg(angle(V_s)-tht_r)
end

%-------------------------------------------------------------------------------------
if select == 8
    % 2 loads, 2 transmission lines
    z_line = 0.5 + 1j*0.3;
    z_1 = 8 + 1j*12;
    z_2 = 6 + 1j*3;
    
    % seen by source, Z
    z_all = z_line + f_para(z_line + z_2, z_1);
    f_mdri("Z_src", z_all, 1);
    pf_all = cos(angle(z_all));
    fprintf("pf = %0.2f  ,  lag by pos angle\n", pf_all);
    
    % reactance needed to bring pf to .95 lag
    pf_new = 0.95;
    tht_new = acos(pf_new);
    
    syms xc;
    eqn = angle(f_para(z_all, xc));
    Xc = double(solve(eqn == tht_new, xc));
    fprintf("need to add reactance:  %f\n", Xc);
    S = f_para(z_all, Xc);
    check = double(angle(S) - tht_new)  % should be 0
end

%-------------------------------------------------------------------------------------
if select == 9
    syms R;
    syms C;
    syms w;
    syms L;
    Z = f_para(R+1j*w*L, 1/(1j*w*C));
    Z_imag = imag(Z); % imag part must be 0...for unity
    eqn = solve(Z_imag == 0, C)
end

%-------------------------------------------------------------------------------------
if select == 10
    v_s = 100;  % 0 deg
    
    z1_mag = 100; % VA
    pf_1 = 0.6; % lag
    z1_r = z1_mag * pf_1;
    tht1 = acos(pf_1);
    z1_i = z1_mag * sin(tht1);
    z1 = z1_r + 1j*z1_i; 
    i_1 = v_s / z1;
    s_1 = v_s * conj(i_1);
    f_mdri("S1", s_1, 1);
    fprintf("|i1|= %0.1f  ,  r1=  %0.1f  ,  X1=  %0.1f\n\n", abs(i_1), real(z1), imag(z1));
    
    
    z2_mag = 70; % VA
    pf_2 = 0.75; %lead
    z2_r = z2_mag * pf_2;
    tht2 = -1 * acos(pf_2);
    z2_i = z2_mag * sin(tht2);
    z2 = z2_r + 1j*z2_i;
    i_2 = v_s / z2;
    s_2 = v_s * conj(i_2);
    f_mdri("S2", s_2, 1);
    fprintf("|i2|= %0.1f  ,  r2=  %0.1f  ,  X2=  %0.1f\n\n", abs(i_2), real(z2), imag(z2));
    
    z3_r = 45; % W
    pf_3 = 0.95; % lag
    tht3 = acos(pf_3);
    z3_mag = z3_r / pf_3;
    z3_i = z3_mag * sin(tht3);
    z3 = z3_r + 1j*z3_i;
    i_3 = v_s / z3;
    s_3 = v_s * conj(i_3);
    f_mdri("S3", s_3, 1);
    fprintf("|i3|= %0.1f  ,  r3=  %0.1f  ,  X3=  %0.1f\n\n", abs(i_3), real(z3), imag(z3));
    
    z_all = f_para(z3, f_para(z1, z2));
    f_mdri("Z", z_all, 1);
    Ii = v_s / z_all;
    check = Ii - i_1 - i_2 - i_3% book is wrong
end

%-------------------------------------------------------------------------------------
if select == 11
    
end

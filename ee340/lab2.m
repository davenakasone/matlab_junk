%{
    1  :  q1
    2  :  q2
    3  :  qb1
    4  :  c1
%}
clear all;
clc;
select = 5;

%-------------------------------------------------------------------------------------
if select == 1
    v_a = 120*exp(1j*deg2rad(0));
    v_b = 120*exp(1j*deg2rad(-120));
    v_c = 120*exp(1j*deg2rad(120));
    z_r = 600*exp(1j*deg2rad(0));
    z_c = 600*exp(1j*deg2rad(-90));
    
    syms v_n;
    eqn = ((v_a-v_n)/z_c) + ((v_b-v_n)/z_r) + ((v_c-v_n)/z_r);
    Vn = double(solve(eqn==0, v_n));
    fprintf("\n");
    f_mdri("Vn", Vn, 1);
    
    Van = v_a - Vn;
    fprintf("\n");
    f_mdri("V_AN", Van, 1);
    Vbn = v_b - Vn;
    fprintf("\n");
    f_mdri("V_BN", Vbn, 1);
    Vcn = v_c - Vn;
    fprintf("\n");
    f_mdri("V_CN", Vcn, 1);
end

%-------------------------------------------------------------------------------------
if select == 2
    v_a = 120*exp(1j*deg2rad(0));
    v_b = 120*exp(1j*deg2rad(120));
    v_c = 120*exp(1j*deg2rad(-120));
    z_r = 600*exp(1j*deg2rad(0));
    z_c = 600*exp(1j*deg2rad(-90));
    
    syms v_n;
    eqn = ((v_a-v_n)/z_c) + ((v_b-v_n)/z_r) + ((v_c-v_n)/z_r);
    Vn = double(solve(eqn==0, v_n));
    fprintf("\n");
    f_mdri("Vn", Vn, 1);
    
    Van = v_a - Vn;
    fprintf("\n");
    f_mdri("V_AN", Van, 1);
    Vbn = v_b - Vn;
    fprintf("\n");
    f_mdri("V_BN", Vbn, 1);
    Vcn = v_c - Vn;
    fprintf("\n");
    f_mdri("V_CN", Vcn, 1);
end

%-------------------------------------------------------------------------------------
if select == 3
    vv = 120;
    v1 = vv * exp(1j*deg2rad(0));
    v2 = vv * exp(1j*deg2rad(-120));
    v3 = vv * exp(1j*deg2rad(120));
    r1 = 170;
    r2 = 240;
    r3 = 300;
    v12 = v1 - v2;
    v23 = v2 - v3;
    v31 = v3 - v1;
    i12 = v12 / r1;
    i23 = v23 / r2;
    i31 = v31 / r3;
    
    s1 = v12 * conj(i12);
    s2 = v23 * conj(i23);
    s3 = v31 * conj(i31);
    S = s1 + s2 + s3;
    fprintf("\n\tpower to load phases:\n\n");
    f_mdri("S1", s1, 1);
    f_mdri("S2", s2, 1);
    f_mdri("S3", s3, 1);
    f_mdri("S", S, 1);
    fprintf("\n\tconfirming by wattmeter, the total power is the same:\n\n");
    
    I1 = i12 - i31;
    I2 = i23 - i12;
    I3 = i31 - i23;
    s12 = (v1-v2) * conj(I1);
    s32 = (v3-v2) * conj(I3);
    SS = s12 + s32;
    f_mdri("s12", s12, 1);
    f_mdri("s32", s32, 1);
    f_mdri("S", SS, 1);
    
    s1 = v1 * conj(I1);
    s2 = v2 * conj(I2);
    s3 = v3 * conj(I3);
    S = s1 + s2 + s3;
    fprintf("\n\tpower supplied by phases:\n\n");
    f_mdri("S1", s1, 1);
    f_mdri("S2", s2, 1);
    f_mdri("S3", s3, 1);
    f_mdri("S", S, 1);
end

%-------------------------------------------------------------------------------------
if select == 4
    v_a = 120 * exp(1j*deg2rad(0));
    v_b = 120 * exp(1j*deg2rad(-120));
    v_c = 120 * exp(1j*deg2rad(120));
    R = 100;
    z_L = 1j*R*sqrt(3);
    z_C = -1j*R*sqrt(3);
    v_ab = v_a - v_b;
    v_bc = v_b - v_c;
    v_ca = v_c - v_a;
    i_12 = v_ab / R; 
    i_23 = v_bc / R; % no inductor
    i_31 = v_ca / R; % no capacitor
    i_a = i_12 - i_31;
    i_b = i_23 - i_12;
    i_c = i_31 - i_23;
    s_a = v_a * conj(i_a);
    s_b = v_b * conj(i_b);
    s_c = v_c * conj(i_c);
    S = s_a + s_b + s_c;
    f_mdri("sa", s_a, 1);
    f_mdri("sb", s_b, 1);
    f_mdri("sc", s_c, 1);
    f_mdri("S", S, 1);
end

%-------------------------------------------------------------------------------------
if select == 3
    vv = 120;
    v1 = vv * exp(1j*deg2rad(0));
    v2 = vv * exp(1j*deg2rad(-120));
    v3 = vv * exp(1j*deg2rad(120));
    r1 = 170;
    r2 = 240;
    r3 = 300;
    v12 = v1 - v2;
    v23 = v2 - v3;
    v31 = v3 - v1;
    i12 = v12 / r1;
    i23 = v23 / r2;
    i31 = v31 / r3;
    
    s1 = v12 * conj(i12);
    s2 = v23 * conj(i23);
    s3 = v31 * conj(i31);
    S = s1 + s2 + s3;
    fprintf("\n\tpower to load phases:\n\n");
    f_mdri("S1", s1, 1);
    f_mdri("S2", s2, 1);
    f_mdri("S3", s3, 1);
    f_mdri("S", S, 1);
    fprintf("\n\tconfirming by wattmeter, the total power is the same:\n\n");
    
    I1 = i12 - i31;
    I2 = i23 - i12;
    I3 = i31 - i23;
    s12 = (v1-v2) * conj(I1);
    s32 = (v3-v2) * conj(I3);
    SS = s12 + s32;
    f_mdri("s12", s12, 1);
    f_mdri("s32", s32, 1);
    f_mdri("S", SS, 1);
    
    s1 = v1 * conj(I1);
    s2 = v2 * conj(I2);
    s3 = v3 * conj(I3);
    S = s1 + s2 + s3;
    fprintf("\n\tpower supplied by phases:\n\n");
    f_mdri("S1", s1, 1);
    f_mdri("S2", s2, 1);
    f_mdri("S3", s3, 1);
    f_mdri("S", S, 1);
end

%-------------------------------------------------------------------------------------
if select == 4
    v_a = 120 * exp(1j*deg2rad(0));
    v_b = 120 * exp(1j*deg2rad(-120));
    v_c = 120 * exp(1j*deg2rad(120));
    R = 100;
    z_L = 1j*R*sqrt(3);
    z_C = -1j*R*sqrt(3);
    v_ab = v_a - v_b;
    v_bc = v_b - v_c;
    v_ca = v_c - v_a;
    i_12 = v_ab / R; 
    i_23 = v_bc / R; % no inductor
    i_31 = v_ca / R; % no capacitor
    i_a = i_12 - i_31;
    i_b = i_23 - i_12;
    i_c = i_31 - i_23;
    s_a = v_a * conj(i_a);
    s_b = v_b * conj(i_b);
    s_c = v_c * conj(i_c);
    S = s_a + s_b + s_c;
    f_mdri("sa", s_a, 1);
    f_mdri("sb", s_b, 1);
    f_mdri("sc", s_c, 1);
    f_mdri("S", S, 1);
end

%-------------------------------------------------------------------------------------
if select == 5
    v_a = 120 * exp(1j*deg2rad(0));
    v_b = 120 * exp(1j*deg2rad(-120));
    v_c = 120 * exp(1j*deg2rad(120));
    R = 100;
    z_L = 1j*R*sqrt(3);
    z_C = -1j*R*sqrt(3);
    v_ab = v_a - v_b;
    v_bc = v_b - v_c;
    v_ca = v_c - v_a;
    i_12 = v_ab / R; 
    i_23 = v_bc / z_C; % capacitor
    i_31 = v_ca / z_L; % inductor
    i_a = i_12 - i_31;
    i_b = i_23 - i_12;
    i_c = i_31 - i_23;
    s_a = v_a * conj(i_a);
    s_b = v_b * conj(i_b);
    s_c = v_c * conj(i_c);
    S = s_a + s_b + s_c;
    f_mdri("sa", s_a, 1);
    f_mdri("sb", s_b, 1);
    f_mdri("sc", s_c, 1);
    f_mdri("S", S, 1);
end
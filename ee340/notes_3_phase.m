%{
    1  : ex1, 3-phase pow measure, 3 watt meters    imbalance is % of max / avg
    2  : ex2, 3-phase pow measure, 2 watt meters
    3  : ex3, 3-phase pow measure, 2 watt meters, delta load
    4  : check f_delta2y and f_y2delta
    5  : ex_balanced, convert delta-->Y
    6  : conversions on source and impedances...
    7  : test f_hp2watts() and f_watts2hp()
    8  : pp, p25
    9  : pp, p26, easy kvl
%}
clc;
select = 9;

%-------------------------------------------------------------------------------------
if select == 1
    v_s = 120;
    r_1 = 100;
    r_2 = 200;
    r_3 = 300;
    
    A = [ r_1, -r_2, 0;
          -r_2, r_2+r_3, -r_3;
          0, -r_3, r_3];
    Y = [0;0;v_s];
    I = A \ Y;
    i_1 = I(1,1);
    i_2 = I(2,1) - I(1,1);
    i_3 = I(3,1) - I(2,1);
    fprintf("\n\ti1: %0.1f  ,  i2: %0.1f  ,  i3: %0.1f\n", i_1, i_2, i_3);
    pow_1 = i_1^2 * r_1;
    pow_2 = i_2^2 * r_2;
    pow_3 = i_3^2 * r_3;
    fprintf("\n\tpow1:  %0.1f  ,  pow2:  %0.2f  ,  pow3:  %0.2f\n", pow_1, pow_2, pow_3);
    i_n = I(3,1);
    fprintf("\n\tIn = %0.1f\n", i_n);
end

%-------------------------------------------------------------------------------------
if select == 2
    v_s = 120;
    r_1 = 100;
    r_2 = 200;
    r_3 = 300;
    
    vp_1 = v_s;
    f_mdri("vp1", vp_1, 1);
    vp_2 = v_s * exp(1j*deg2rad(-120));
    f_mdri("vp2", vp_2, 1);
    vp_3 = v_s * exp(1j*deg2rad(120));
    f_mdri("vp3", vp_3, 1);
    fprintf("\n");
    
    vl_12 = vp_1 - vp_2;
    f_mdri("vl_12", vl_12, 1);
    vl_32 = vp_3 - vp_2;
    f_mdri("vl_32", vl_32, 1);
    v_n = vl_12 - vl_32;
    f_mdri("Vn", v_n, 1);
    fprintf("\n");
    
    % matrix don't really work...
    A = [r_1+r_2, -r_2; -r_2, r_2 + r_3];
    Y = [vp_1 - vp_2; vp_2 - vp_3];
    I = A \ Y;
    i_1 = I(1,1);
    f_mdri("i1", i_1, 1);
    i_2 = I(1,1) - I(2,1);
    f_mdri("i2", i_2, 1);
    i_3 = -1 * I(2,1);
    f_mdri("i3", i_3, 1);
    fprintf("\n");
    
    % KCL on neutral is better
    vn = ((1/r_1)+(1/r_2)+(1/r_3))^(-1) * ((vp_1/r_1)+(vp_2/r_2)+(vp_3/r_3));
    i_1 = (vp_1-vn)/r_1;
    f_mdri("i1", i_1, 1);
    i_2 = (vp_2-vn)/r_2;
    f_mdri("i2", i_2, 1);
    i_3 = (vp_3-vn)/r_3;
    f_mdri("i3", i_3, 1);
    fprintf("\n");
    
    pow_12 = vl_12 * conj(i_1);
    f_mdri("pow12", pow_12, 1);
    pow_32 = vl_32 * conj(i_3);
    f_mdri("pow32", pow_32, 1);
    pow_tot =  pow_12 + pow_32;
    f_mdri("pow_tot", pow_tot, 1);
    fprintf("\n");
    pow_1 = vp_1 * conj(i_1);
    f_mdri("pow1", pow_1, 1);
    pow_2 = vp_2 * conj(i_2);
    f_mdri("pow2", pow_2, 1);
    pow_3 = vp_3 * conj(i_3);
    f_mdri("pow3", pow_3, 1);
    f_mdri("pow_tot", pow_1+pow_2+pow_3, 1);
    fprintf("\n");
    
    i_avg = (abs(i_1)+abs(i_2)+abs(i_3))/3;
    i_diff = abs(abs([i_1, i_2, i_3]) - i_avg);
    i_imb = max(i_diff, [], 'all') / i_avg;
    fprintf("\tmax I_imbalance:  %0.1f percent\n", i_imb*100);
end

%-------------------------------------------------------------------------------------
if select == 3
    v_s = 120;
    r_1 = 100;
    r_2 = 200;
    r_3 = 300;
    
    vp_1 = v_s;
    f_mdri("vp1", vp_1, 1);
    vp_2 = v_s * exp(1j*deg2rad(-120));
    f_mdri("vp2", vp_2, 1);
    vp_3 = v_s * exp(1j*deg2rad(120));
    f_mdri("vp3", vp_3, 1);
    fprintf("\n");
    
    vl_12 = vp_1 - vp_2;
    f_mdri("vl_12", vl_12, 1);
    vl_32 = vp_3 - vp_2;
    f_mdri("vl_32", vl_32, 1);
    v_n = vl_12 - vl_32;
    f_mdri("Vn", v_n, 1);
    fprintf("\n");
    
    A = [r_1, 0, -r_1; 0, r_2, -r_2; -r_1, -r_2, r_1+r_2+r_3];
    Y = [vp_1-vp_2; vp_2-vp_3; 0];
    I = A\Y;
    i_12 = I(1,1)-I(3,1);
    f_mdri("i12", i_12, 1);
    i_23 = I(2,1)-I(3,1);
    f_mdri("i23", i_23, 1);
    i_31 = -I(3,1);
    f_mdri("i31", i_31, 1);
    fprintf("\n");
    
    i_1 = I(1,1);
    f_mdri("i1", i_1, 1);
    i_3 = -I(2,1);
    f_mdri("i2", i_3, 1);
    fprintf("\n");
    
    pow_12 = vl_12 * conj(i_1);
    f_mdri("pow_12", pow_12, 1);
    pow_32 = vl_32 * conj(i_3);
    f_mdri("pow_32", pow_32, 1);
    pow_tot = pow_12 + pow_32;
    f_mdri("pow_tot", pow_tot, 1);
end

%-------------------------------------------------------------------------------------
if select == 4
    dela = 10;
    delb = 20;
    delc = 30;
    [y1, y2, y3] = f_delta2y(dela, delb, delc);
    fprintf("dela: %0.2f  ,  delb: %0.2f  ,  delc: %0.2f\n", dela, delb, delc);
    fprintf("y1: %0.2f  ,  y2: %0.2f  ,  y3: %0.2f\n", y1, y2, y3);
    [dela, delb, delc] = f_y2delta(y1, y2, y3);
    fprintf("dela: %0.2f  ,  delb: %0.2f  ,  delc: %0.2f\n", dela, delb, delc);
end

%-------------------------------------------------------------------------------------
if select == 5
    vl_1 = 170;
    vl_2 = 170 * exp(1j*deg2rad(-120));
    vl_3 = 170 * exp(1j*deg2rad(120));
    z_p = 20 + 1j*5; % balanced, all deltas have this Z
    [any, any1, any2] = f_delta2y(z_p, z_p, z_p); % should all be the same
    Z = any;
    i_1 = vl_1 / Z;
end

%-------------------------------------------------------------------------------------
if select == 6
    z_L = (1/10)*1j;
    z_c = 1/(1j*10);
    [zcy, zcyy, zcyyy] = f_delta2y(z_c, z_c, z_c); % all the same
    vs_la = 1;
    vs_lb = exp(1j*deg2rad(-120));
    vs_lc = exp(1j*deg2rad(120));
    %vs_dra = 1;
    %vs_drb = exp(1j*deg2rad(-120));
    %vs_drc = exp(1j*deg2rad(120));
    %[vs_r1, vs_r2, vs_r3] = f_delta2y(vs_dra, vs_drb, vs_drc) bad idea
    [vs_r1, vs_r2, vs_r3] = f_delta2y(1, 1, 1); % all the same, get angles on them
    vs_yr1 = vs_r1 * exp(1j*deg2rad(0));
    vs_yr2 = vs_r2 * exp(1j*deg2rad(-120));
    vs_yr3 = vs_r3 * exp(1j*deg2rad(120));
    
    node_a = ( (vs_la/z_L) + (vs_yr1/z_L) ) / ( (2/z_L) + (1/zcy))
end

%-------------------------------------------------------------------------------------
if select == 7
    h_p = 100;
    h2w = f_hp2watts(h_p)
    wat = f_watts2hp(h2w)
end

%-------------------------------------------------------------------------------------
if select == 8
    v_a = 120;
    v_b = 120 * exp(1j*deg2rad(-120));
    v_c = 120 * exp(1j*deg2rad(120));
    z_ab = 60;
    z_bc = 90;
    
    i_1 = (v_a - v_b) / z_ab;
    f_mdri("i1", i_1, 1);
    p_1 = (v_a - v_b) * conj(i_1);
    f_mdri("pow1", p_1, 1);
    i_2 = (v_c - v_b) / z_bc;
    f_mdri("i2", i_2, 1);
    p_2 = (v_c - v_b) * conj(i_2);
    f_mdri("pow2", p_2, 1);
    i_b = -(i_1 + i_2);
    fprintf("\tmiddle line, ib =  %0.2f  A\n\n", abs(i_b));
    
    % now r_ca is not open
    z_ac = 120;
    A = [z_ab, 0, -z_ab; 0, z_bc, -z_bc; -z_ab, -z_bc, z_ab+z_bc+z_ac];
    Y = [v_a-v_b; v_b-v_c; 0];
    I = A \ Y;
    i_a = I(1,1);
    i_b = I(2,1) - I(1,1);
    i_c = -I(2,1);
    f_mdri("ia", i_a, 1);
    p_a = (v_a - v_b) * conj(i_a);
    f_mdri("p1", p_a, 1);
    f_mdri("ic", i_c, 1);
    p_b = (v_c-v_b) * conj(i_c);
    f_mdri("p2", p_b, 1);
    fprintf("\tmiddle line, ib =  %0.2f  A\n\n", abs(i_b));
end

%-------------------------------------------------------------------------------------
if select == 9
    v_a = 120 * exp(1j*deg2rad(0));
    v_b = 120 * exp(1j*deg2rad(-120));
    v_c = 120 * exp(1j*deg2rad(120));
    R = 3;
    z_ab = 10;
    z_bc = 5;
    z_ca = 10-1j*5;
    A = [R+z_ab+R, -R, -z_ab; -R, R+z_bc+R, -z_bc; -z_ab, -z_bc, z_ab+z_bc+z_ca];
    Y = [v_a-v_b; v_b-v_c; 0];
    I = A \ Y;
    i_1 = I(1,1);
    i_2 = I(1,1) - I(2,1);
    i_3 = -I(2,1);
    f_mdri("i1", i_1, 1);
    f_mdri("i2", i_2, 1);
    f_mdri("i3", i_3, 1);
end
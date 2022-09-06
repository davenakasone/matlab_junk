%{
    1: p1, solve powers provided by phase, no balancing yet
    2: p2, balance the lo
    3: p2, again
%}
close all;
clear;
clc;
select = 1;

%-------------------------------------------------------------------------------------
if select == 1
    V_a = 120 * exp(1j*deg2rad(0));
    V_b = 120 * exp(1j*deg2rad(-120));
    V_c = 120 * exp(1j*deg2rad(120));
    R_12 = 10;
    R_23 = 15;
    R_31 = 20;
    
    % solve the KVL mesh
    A_mtx = [R_12 , 0   , -R_12; 
             0    , R_23, -R_23;
             -R_12, -R_23, R_12 + R_23 + R_31];
    Y_mtx = [V_a - V_b;
             V_b - V_c;
             0];
    I_mtx = A_mtx \ Y_mtx;
    I_a = I_mtx(1,1);
    I_b = I_mtx(2,1) - I_mtx(1,1);
    I_c = -1 * I_mtx(2,1);
    S_a = V_a * conj(I_a);
    S_b = V_b * conj(I_b);
    S_c = V_c * conj(I_c);
    S_total = S_a + S_b + S_c;
    f_mdri("Sa", S_a, 1);
    f_mdri("Sb", S_b, 1);
    f_mdri("Sc", S_c, 1);
    f_mdri("S_total", S_total, 1);
    
    % check with KCL
    fprintf("\n");
    i_12 = (V_a - V_b) / R_12;
    i_23 = (V_b - V_c) / R_23;
    i_31 = (V_c - V_a) / R_31;
    I_aa = i_12 - i_31;
    I_bb = i_23 - i_12;
    I_cc = i_31 - i_23;
    S_aa = V_a * conj(I_aa);
    S_bb = V_b * conj(I_bb);
    S_cc = V_c * conj(I_cc);
    S_total = S_aa + S_bb + S_cc;
    f_mdri("Sa", S_aa, 1);
    f_mdri("Sb", S_bb, 1);
    f_mdri("Sc", S_cc, 1);
    f_mdri("S_total", S_total, 1);
end


%-------------------------------------------------------------------------------------
if select == 2
    V_a = 120 * exp(1j*deg2rad(0));
    V_b = 120 * exp(1j*deg2rad(-120));
    V_c = 120 * exp(1j*deg2rad(120));
    R_12 = 10;
    R_23 = 15;
    R_31 = 20;
    
    % use the super-position reactive bank
    Z_12 = 1j*sqrt(3) * (R_23 - R_31);
    Z_23 = 1j*sqrt(3) * (R_31 - R_12);
    Z_31 = 1j*sqrt(3) * (R_12 - R_23);
    f_mdri("Z_12", Z_12, 1);
    f_mdri("Z_23", Z_23, 1);
    f_mdri("Z_31", Z_31, 1);
    fprintf("\n");
    
    % get the impedance in parallel
    imp_12 = f_para(R_12, Z_12);
    imp_23 = f_para(R_23, Z_23);
    imp_31 = f_para(R_31, Z_31);
    
    % solve the KVL mesh
    A_mtx = [imp_12 , 0   , -imp_12; 
             0    , imp_23, -imp_23;
             -imp_12, -imp_23, imp_12 + imp_23 + imp_31];
    Y_mtx = [V_a - V_b;
             V_b - V_c;
             0];
    I_mtx = A_mtx \ Y_mtx;
    I_a = I_mtx(1,1);
    I_b = I_mtx(2,1) - I_mtx(1,1);
    I_c = -1 * I_mtx(2,1);
    S_a = V_a * conj(I_a);
    S_b = V_b * conj(I_b);
    S_c = V_c * conj(I_c);
    S_total = S_a + S_b + S_c;
    f_mdri("Sa", S_a, 1);
    f_mdri("Sb", S_b, 1);
    f_mdri("Sc", S_c, 1);
    f_mdri("S_total", S_total, 1);
    
    % check with KCL
    fprintf("\n");
    i_12 = (V_a - V_b) / imp_12;
    i_23 = (V_b - V_c) / imp_23;
    i_31 = (V_c - V_a) / imp_31;
    I_aa = i_12 - i_31;
    I_bb = i_23 - i_12;
    I_cc = i_31 - i_23;
    S_aa = V_a * conj(I_aa);
    S_bb = V_b * conj(I_bb);
    S_cc = V_c * conj(I_cc);
    S_total = S_aa + S_bb + S_cc;
    f_mdri("Sa", S_aa, 1);
    f_mdri("Sb", S_bb, 1);
    f_mdri("Sc", S_cc, 1);
    f_mdri("S_total", S_total, 1);
end


%-------------------------------------------------------------------------------------
if select == 3
    V_a = 120 * exp(1j*deg2rad(0));
    V_b = 120 * exp(1j*deg2rad(-120));
    V_c = 120 * exp(1j*deg2rad(120));
    R_12 = 10;
    R_23 = 15;
    R_31 = 20;
    
    syms z12;
    syms z23;
    syms z31;
    iimp_12 = f_para(R_12, 1j*z12);
    iimp_23 = f_para(R_23, 1j*z23);
    iimp_31 = f_para(R_31, 1j*z31);
    
    A_mtx = [iimp_12 , 0   , -iimp_12; 
             0    , iimp_23, -iimp_23;
             -iimp_12, -iimp_23, iimp_12 + iimp_23 + iimp_31];
    Y_mtx = [V_a - V_b;
             V_b - V_c;
             0];
    I_mtx = A_mtx \ Y_mtx;
    I_a = I_mtx(1,1);
    I_b = I_mtx(2,1) - I_mtx(1,1);
    I_c = -1 * I_mtx(2,1);
    S_a = V_a * conj(I_a);
    S_b = V_b * conj(I_b);
    S_c = V_c * conj(I_c);
    [z_12, z_23, z_31] = vpasolve([imag(S_a)==0, imag(S_b)==0, imag(S_c)==0],...
        [z12, z23, z31]);
    Z_12 = z_12 * 1j
    Z_23 = z_23 * 1j
    Z_31 = z_31 * 1j
    
    % take the new reactive values and put them in parallel
    imp_12 = f_para(R_12, Z_12);
    imp_23 = f_para(R_23, Z_23);
    imp_31 = f_para(R_31, Z_31);
    
    % solve the KVL mesh
    A_mtx = [imp_12 , 0   , -imp_12; 
             0    , imp_23, -imp_23;
             -imp_12, -imp_23, imp_12 + imp_23 + imp_31];
    Y_mtx = [V_a - V_b;
             V_b - V_c;
             0];
    I_mtx = A_mtx \ Y_mtx;
    I_a = I_mtx(1,1);
    I_b = I_mtx(2,1) - I_mtx(1,1);
    I_c = -1 * I_mtx(2,1);
    S_a = V_a * conj(I_a);
    S_b = V_b * conj(I_b);
    S_c = V_c * conj(I_c);
    S_total = S_a + S_b + S_c;
    f_mdri("Sa", S_a, 1);
    f_mdri("Sb", S_b, 1);
    f_mdri("Sc", S_c, 1);
    f_mdri("S_total", S_total, 1);
    
    % check with KCL
    fprintf("\n");
    i_12 = (V_a - V_b) / imp_12;
    i_23 = (V_b - V_c) / imp_23;
    i_31 = (V_c - V_a) / imp_31;
    I_aa = i_12 - i_31;
    I_bb = i_23 - i_12;
    I_cc = i_31 - i_23;
    S_aa = V_a * conj(I_aa);
    S_bb = V_b * conj(I_bb);
    S_cc = V_c * conj(I_cc);
    S_total = S_aa + S_bb + S_cc;
    f_mdri("Sa", S_aa, 1);
    f_mdri("Sb", S_bb, 1);
    f_mdri("Sc", S_cc, 1);
    f_mdri("S_total", S_total, 1);
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END> 
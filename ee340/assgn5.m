%{
    1  : p1.1, island, pow
    2  : p1.2, island, kvl
    3  : p1.3, island, update
    4  : p2.1, balance para
    5  : p2.2, balance para, equal loads
    6  : p2.3, balance para, equal loads, then a change  ...solve()
    7  : p2.4, solve the new set points for #6
    8  : quiz...fuck


    ...infinite bus stays fixed
%}
format compact;
close all;
clear all;
clc;
select = 8;


%-------------------------------------------------------------------------------------
if select == 1
    V_phi = 120; % as rated ...reference to 0
    X_s = 0.5; % as given...  no R_A
    R_load = 4;
    
    I_A = V_phi/R_load;
    fprintf("\nI_A = %0.2f A\n", I_A);
    P = I_A^2 * R_load;
    fprintf("P, per phase = %0.2f W\n", P);
    E_A = V_phi + I_A * 1j*X_s;
    f_mdri("E_A", E_A, 1);
end


%-------------------------------------------------------------------------------------
if select == 2
    V_phi = 120; % as rated ...reference to 0
    X_s = 0.5; % as given...  no R_A
    R_load = f_para(4, 12);
    
    I_A = V_phi/R_load;
    fprintf("\nI_A = %0.2f A\n", I_A);
    P = I_A^2 * R_load;
    fprintf("P, per phase = %0.2f W\n", P);
    E_A = V_phi + I_A * 1j*X_s;
    f_mdri("E_A", E_A, 1);
end


%-------------------------------------------------------------------------------------
if select == 3
    f = 60;
    X_s1 = 0.5;
    E_A1 = 120 + 30*1j*X_s1;
    P_1 = 3600; % from problem 1.1
    R_2 = 3; % from problem 1.2
    
    V_T = sqrt(P_1 * R_2);
    I_A = V_T / R_2;
    fprintf("\nV_T =  %0.2f V\n", V_T);
    
    temp = (P_1 * X_s1) / (V_T * abs(E_A1));
    delta = asin(temp);
    fprintf("delta = % 0.2f\n", rad2deg(delta));
    
    X_s = V_T * tan(delta) / I_A;
    f_mdri("Xs", X_s, 1);
    
    E_A = (P_1 * X_s) / (V_T * sin(delta));
    fprintf("E_A = % 0.2f V\n", E_A);
    
    ff = f * (X_s / X_s1);
    fprintf("f is now  %0.2f Hz\n", ff);
    
    n_m = ff * 120 / 4;
    fprintf("new rpms, 4 pole machine,  %0.2f\n", n_m);
end


%-------------------------------------------------------------------------------------
if select == 4
    f = 60;
    p_1 = 15;
    p_2 = 10;
    p_3 = 5;
    r_1 = 0.1;
    r_2 = 0.15;
    r_3 = 0.2;
    
    b_1 = f + r_1*p_1
    b_2 = f + r_2*p_2
    b_3 = f + r_3*p_3
    
end


%-------------------------------------------------------------------------------------
if select == 5
    f = 60;
    p_1 = 10;
    p_2 = 10;
    p_3 = 10;
    r_1 = 0.1;
    r_2 = 0.15;
    r_3 = 0.2;
    
    b_1 = f + r_1*p_1
    b_2 = f + r_2*p_2
    b_3 = f + r_3*p_3
end


%-------------------------------------------------------------------------------------
if select == 6
    b_1 = 61;   % from 2.2
    b_2 = 61.5; % from 2.2;
    b_3 = 62;   % from 2.2
    r_1 = 0.1;
    r_2 = 0.15;
    r_3 = 0.2;
    new_pow = 45;
    
    syms ff;
    syms p1;
    syms p2;
    syms p3;
    
    eqn_a = r_1 * p1;  % f1 = b1 - r1 * p1 -->    b1 - f1 == r1 * p1
    eqn_aa = b_1 - ff;
    eqn_b = r_2 * p2;
    eqn_bb = b_2 - ff;
    eqn_c = r_3 * p3;
    eqn_cc = b_3 - ff;
    eqn_p = p1 + p2 + p3;
    [f_f, p_1, p_2, p_3] = ...
        solve([eqn_a==eqn_aa, eqn_b==eqn_bb, eqn_c==eqn_cc, eqn_p==new_pow], [ff, p1, p2, p3]);
    fprintf("freq:  %0.2f Hz\n", f_f);
    fprintf("p1:  %0.2f MW\n", p_1);
    fprintf("p2:  %0.2f MW\n", p_2);
    fprintf("p3:  %0.2f MW\n", p_3);
    % alternative     [R] [p] = [b-f], find f, then ps
    A = [r_1, 0, 0; 
         0  , r_2, 0;
         0, 0, r_3];
    Y = [b_1-ff; b_2-ff; b_3-ff];
    pp = A \ Y;
    FF = solve(pp(1)+pp(2)+pp(3)==new_pow, ff);
    fprintf("freq:  %0.2f Hz\n", FF);
    fprintf("p1:  %0.2f MW\n", subs(pp(1), ff, FF));
    fprintf("p2:  %0.2f MW\n", subs(pp(2), ff, FF));
    fprintf("p3:  %0.2f MW\n", subs(pp(3), ff, FF));
end

%-------------------------------------------------------------------------------------
if select == 7
    p_1 = 15;
    p_2 = p_1;
    p_3 = p_2;
    f = 60;
    r_1 = 0.1;
    r_2 = 0.15;
    r_3 = 0.2;
    
    b_1 = f + r_1 * p_1
    b_2 = f + r_2 * p_2
    b_3 = f + r_3 * p_3
end


%-------------------------------------------------------------------------------------
if select == 8
    S_rate = 300e6;
    pf = 0.8;
    f = 60;
    v_ll = 15e3;
    z_in = 0 + 1j*0.2;
    P = 240e6;
    Q = 180e6;
    Pf = 80e6;
    Qf = 60e6;
    Vf = 8.66e3;
    X_s = real(z_in);
    tht = acos(pf); % fuck
    Sf = Pf + 1j*Qf;
    
    I_A = (abs(Sf) / Vf) * pf
    E_A = Vf + abs(I_A)*pf*1j*X_s
    f_mdri("EA", E_A, 1);
    
    
end


%-------------------------------------------------------------------------------------
if select == 99
    
end
%%%%%%%%~~~~~~~~~END> 
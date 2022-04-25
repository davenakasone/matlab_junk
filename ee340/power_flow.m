%{
    power flow
    1  :  ex, slide 5
    2  :  ex, slide, limits
    3  :  ex, quiz
%}
format compact;
close all;
clear;
clc;
select = 3;

%-------------------------------------------------------------------------------------
if select == 1
    % make a per-phase equivelant circuit and use KCL
    % generator 1, transformer 1, a line, transformer 2, a motor
    S_base = 100e6; % of G1
    V_base = 13.8e3;
    
    g1_rate_S = 13.8e3;
    g1_pu_R = 0.1;
    g1_pu_Xs = 0.05;
    
    t1_rate_S = 100e6;
    t1_rate_V = (13.8/110)*10^3;
    t1_pu_R = 0.01;
    t1_pi_Xs = 0.05;

    t2_rate_S = 50e6;
    t2_rate_V = (120/14.4)*10^3;
    t2_rate_R = 0.01;
    t2_rate_Xs = 0.05;

    m_rate_S = 50e6;
    m_rate_V = 13.8e3;
    m_pu_R = 0.1;
    m_pu_Xs = 1.1;

    line_R = 15;
    line_X = 75;

    V_1_base = 13.8e3;
    V_2_base = V_1_base * (t1_rate_S/t1_rate_V)
end


%-------------------------------------------------------------------------------------
if select == 2
    p_demand = 975;

    syms p_1;
    alpha_1 = 500;
    beta_1 = 5.3;
    gamma_1 = 0.004;
    p1_ll = 200;
    p1_ul = 450;
    c_1 = alpha_1 + beta_1 * p_1 + gamma_1 * p_1.^2;

    syms p_2;
    alpha_2 = 400;
    beta_2 = 5.5;
    gamma_2 = 0.006;
    p2_ll = 150;
    p2_ul = 350;
    c_2 = alpha_2 + beta_2 * p_2 + gamma_2 * p_2.^2;

    syms p_3;
    alpha_3 = 200;
    beta_3 = 5.8;
    gamma_3 = 0.009;
    p3_ll = 100;
    p3_ul = 225;
    c_3 = alpha_3 + beta_3 * p_3 + gamma_3 * p_3.^2;

    eqn_pow = p_1 + p_2 + p_3
    eqn_lam1 = diff(c_1, p_1)
    eqn_lam2 = diff(c_2, p_2)
    eqn_lam3 = diff(c_3, p_3)

    lamb_den =  (1/(2*gamma_1)) + (1/(2*gamma_2)) + (1/(2*gamma_3)); 
    lamb_num = (beta_1/(2*gamma_1)) + (beta_2/(2*gamma_2)) + (beta_3/(2*gamma_3));
    lamb = (p_demand + lamb_num) / lamb_den % $/MW
    
    P1 = (lamb - beta_1)/(2*gamma_1)
    P2 = (lamb - beta_2)/(2*gamma_2)
    P3 = (lamb - beta_3)/(2*gamma_3)
end


%-------------------------------------------------------------------------------------
if select == 3
    p_demand = 85;

    syms p_1;
    alpha_1 = 200;
    beta_1 = 40;
    gamma_1 = 0.1;
    p1_ll = 20;
    p1_ul = 125;
    c_1 = alpha_1 + beta_1 * p_1 + gamma_1 * p_1.^2;

    syms p_2;
    alpha_2 = 150;
    beta_2 = 30;
    gamma_2 = 0.125;
    p2_ll = 20;
    p2_ul = 125;
    c_2 = alpha_2 + beta_2 * p_2 + gamma_2 * p_2.^2;

    eqn_pow = p_1 + p_2
    eqn_lam1 = diff(c_1, p_1)
    eqn_lam2 = diff(c_2, p_2)

    lamb_den =  (1/(2*gamma_1)) + (1/(2*gamma_2)); 
    lamb_num = (beta_1/(2*gamma_1)) + (beta_2/(2*gamma_2));
    fprintf("\nQ1\n");
    lamb = (p_demand + lamb_num) / lamb_den % $/MW
    P1 = (lamb - beta_1)/(2*gamma_1)
    P2 = (lamb - beta_2)/(2*gamma_2)

    fprintf("\nQ2\n");
    syms p_p;
    c_total = c_1 + c_2;
    c_c = subs(c_total, p_1, p_p);
    c_cc = subs(c_c, p_2, p_p)
    c_ccc = diff(c_cc, p_p);
    lamm = (p_demand + (70/(2*(9/40))))/(1/(2*(9/40)))

    fprintf("\nQ3\n");
    p_11 = (55-beta_1)/(2*gamma_1)
    p_22 = (55-beta_2)/(2*gamma_2)
    ptot = p_11 + p_22
    


end


%-------------------------------------------------------------------------------------
if select == 99
    
end
%%%%%%%%~~~~~~~~~END> 
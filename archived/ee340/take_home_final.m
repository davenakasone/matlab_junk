%{
    take-home-final
    
    1  :  solve the economic dispatch
%}
format compact;
close all;
clear;
clc;
select = 2;


%-------------------------------------------------------------------------------------
if select == 1
    P_load_1 = 70;
    Q_load_1 = 15;
    alpha_1 = 500;
    beta_1 = 20;
    gamma_1 = 0.02;
    P1_ll = 20;
    P1_ul = 200;

    P_load_3 = 110;
    Q_load_3 = 40;
    alpha_3 = 400;
    beta_3 = 22;
    gamma_3 = 0.02;
    P3_ll = 50;
    P3_ul = 150;

    P_load_4 = 90;
    Q_load_4 = 20;
    alpha_4 = 300;
    beta_4 = 24;
    gamma_4 = 0.02;
    P4_ll = 30;
    P4_ul = 120;

    P_load_total = P_load_1 + P_load_3 + P_load_4;
    Q_load_total = Q_load_1 + Q_load_3 + Q_load_4;
    fprintf("\nthe load must equal the demad:\n");
    fprintf("P_total:  %0.1f MW\n", P_load_total);
    fprintf("Q_total:  %0.1f MVAR\n", Q_load_total);

    temp_a = (beta_1/(2*gamma_1)) + (beta_3/(2*gamma_3)) + (beta_4/(2*gamma_4));
    temp_b = (1/(2*gamma_1)) + (1/(2*gamma_3)) + (1/(2*gamma_4));
    lambda = (P_load_total + temp_a) / temp_b;
    fprintf("\nlambda:  %0.2f $/MWhr\n", lambda);

    P_1 = (lambda - beta_1) / (2 * gamma_1);
    P_3 = (lambda - beta_3) / (2 * gamma_3);
    P_4 = (lambda - beta_4) / (2 * gamma_4);
    fprintf("\nP1 [ %0.1d : %0.1d ]  ,  operate=  %0.1f MW\n", P1_ll, P1_ul, P_1);
    fprintf("P3 [ %0.1d : %0.1d ]  ,  operate=  %0.1f MW\n", P3_ll, P3_ul, P_3);
    fprintf("P4 [ %0.1d : %0.1d ]  ,  operate=  %0.1f MW\n", P4_ll, P4_ul, P_4);
    fprintf("check, P_demand:  %0.1f MW ,  P_supply:  %0.1f MW\n",...
        P_load_total, P_1 + P_3 + P_4);
end

%-------------------------------------------------------------------------------------
if select == 2 % quiz8
    eqn = [4, 48, 96, 0];
    rtz = roots(eqn)
    
end
%%%%%%%%~~~~~~~~~END>
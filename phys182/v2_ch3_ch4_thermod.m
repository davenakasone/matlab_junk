%{
    v2, ch3, ch4, thermodynamics laws #1, #2
    
    2  :  p2, ideal gas and isothermic expansion
    3  :  p3, engine and power
    4  :  p4, work, using area in an enclosed contour
    5  :  p5, compare efficiency
    6  :  p6, basic Carnot engine
    7  :  p7, efficiency evaluation
    8  :  p8, gas engine, Otto cycle, monatomic ideal gas (todo)
    9  :  p9, refergerator
    10 :  p10, hvac
    11 :  p11, carnot engine
    12 :  p12, entropy
    13 :  p13, entropy
    14 :  p14, gas
%}
clc;
clear;
select = 8;


%%%%%~~~~
if select == 2
    R = 8.314; % gas constant "R" in [J/(mol K)]
    n = 1; % moles of ideal gas 
    work = 3800; %@@@ J
    % expands isothermally (const T)
    p_final = 1e5; % Pa
    v_final = 0.035; %@@@ m^3

    v_initial = v_final * exp(-1 * work / (p_final * v_final));
    fprintf("a) v_initial = %0.5f\n", v_initial);

    T = (p_final * v_final) / (n * R);
    fprintf("b) T = %0.3f K\n", T);
end


%%%%%~~~~
if select == 3
    Q_hot = 1845; %@@@ J
    Q_cold = 1340; %@@@ J
    t_cycle = 0.21;%@@@ s
    W_actual = 240; %@@@ J

    eta_ideal = 1 - (Q_cold/Q_hot);
    fprintf("a)  eta=  %0.3f\n", eta_ideal);

    W_ideal = Q_hot - Q_cold;
    fprintf("b) W = %0.1f J\n", W_ideal);

    P_ideal = W_ideal / t_cycle;
    fprintf("c) P= %0.2f W\n", P_ideal);

    eta_actual = W_actual / Q_hot;
    fprintf("d) W= %0.4f\n", eta_actual);

    P_actual = W_actual / t_cycle;
    fprintf("e) P=  %0.3f W\n", P_actual);
end


%%%%%~~~~
if select == 4
    V_A = 0.55e-3;
    V_B = 3.8e-3;
    p_A = 2.85e6;
    p_B = 1.86e6;
    p_C = 0.55e6;
    p_D = 1.08e6;

    r_1 = (1/2) * (p_B - p_C) * (V_B - V_A);
    r_2 = (1/2) * (p_A - p_D) * (V_B - V_A);
    W = r_1 + r_2;
    fprintf("W = %0.1f J\n", W);
end


%%%%%~~~~
if select == 5
    eta_norm = 29/100; %@@@
    eta_imp = 43/100; %@@@
    Q_in = 2.75e14; %@@@ J, heat transfer in for 1 day

    W_norm = eta_norm * Q_in;
    W_imp = eta_imp * Q_in;
    W_diff = W_imp - W_norm

    Q_diff = -1 * W_diff
end


%%%%%~~~~
if select == 6
    T_h = 910 + 273; %@@@ --> K
    T_c = 520 + 273; %@@@ --> K
    eta = 1 - (T_c/T_h);
    fprintf("eta=  %0.4f\n", eta);

end


%%%%%~~~~
if select == 7
    T_h = 330 + 273; %@@@ --> K
    T_c = 265 + 273; %@@@ --> K
    eta_claimed = 0.37;

    eta_max = (T_h - T_c) / T_h;
    if eta_max > eta_claimed
        fprintf("claimed eta=  %0.4f  <  eta,max=  %0.4f  ...could work\n", eta_claimed, eta_max);
    else
        fprintf("claimed eta=  %0.4f  >  eta,max=  %0.4f  ...not possible\n", eta_claimed, eta_max);
    end
    
end


%%%%%~~~~
if select == 8
    p_1 = 1.1e5; %@@@ Pa
    V_1 = 0.035; %@@@ m^3
    T_1 = 295; %@@@ K
    aa = 9.5; %@@@   V2 = V1/aa
    bb = 1.56; %@@@  T3 = T2/bb
    V_2 = V_1/aa;


    nR = p_1 * V_1 / T_1;
    fprintf("a) nR=  %0.4f\n", nR);
    fprintf("b) W12 = (3/2)(p1V1((V1/V2)^(2/3) - 1)\n");

    W_12 = (3/2) * p_1 * V_1 * ((V_1/V_2)^(2/3) - 1);
    fprintf("c) W12=  %0.3f J\n", W_12);

    T_2 = T_1 * (aa^(2/3));
    fprintf("d) T2=  %0.3f K\n", T_2);

    Qh = (3/2) * p_1 * V_1 * (aa^(2/3)) * (bb-1);
    fprintf("e)  Qh=  %0.3f\n", Qh);

    p_3 = p_1 * bb * aa^(5/3);
    fprintf("f) p3=  %d\n", p_3);

    W_34 = p_1 * V_1 * bb * (1 - aa^(2/3))/(2/3);
    fprintf("g) W34=  %0.3f\n", W_34);

    Qc = abs((3/2) * p_1 * V_1 * (1-bb));
    fprintf("h) Qc=  %0.1f\n", Qc);

    eta = 100 * (1 - 1/(aa^(2/3)));
    fprintf("i) eta = %0.3f %%\n", eta);

end


%%%%%~~~~
if select == 9
    T_c = -9 + 273; %@@@ --> K
    T_h = 27 + 273; %@@@ --> K
    W = 2.75e3; %@@@ J
    Q_c = 12.4e3; %@@@ J/s
    m = 2.8; %@@@ kg
    T_0 = 3.5 + 273; %@@@ --> K
    T_1 = -1 + 273; % --> K
    c = 4000; % J/(kg K)
    Lf = 280e3; % J/Kg

    eta_max = T_c / (T_h - T_c);
    fprintf("a) Kmax=  %0.3f\n", eta_max);

    eta = Q_c/W;
    fprintf("b) eta=  %0.3f\n", eta);

    Q_reduced1 = m * c * (T_0 - T_1);
    dt1 = Q_reduced1 / Q_c;
    fprintf("c) dt1=  %0.3f s\n", dt1);

    Q_reduced2 = m * Lf;
    dt2 = Q_reduced2 / Q_c;
    fprintf("d) dt2=  %0.3f s\n", dt2);

end


%%%%%~~~~
if select == 10
    T_out = 37 + 273; %@@@ --> K
    T_21 = 21 + 273; % --> K
    T_25 = 25 + 273; % --> K

    p_21 = (T_out - T_21)^2 / T_21;
    p_25 = (T_out - T_25)^2 / T_25;
    A = p_21 / p_25;
    fprintf("rate must be:  %0.3f\n", A);
end


%%%%%~~~~
if select == 11
    eta = 0.45; %@@@
    cop = (1/eta) - 1;
    fprintf("COP=  %0.3f\n", cop);
    fprintf("COP=  %0.3f\n", 1/eta);
end


%%%%%~~~~
if select == 12
    m = 2.1; %@@@ kg
    T_1 = 18.5 + 273.15; %@@@ --> K
    T_2 = T_1 + 10;
    cp = 1000;

    dS = m * cp * log(T_2/T_1);
    fprintf("dS=  %0.3f cal/K\n", dS);
end


%%%%%~~~~
if select == 13
    m = 770/1000; %@@@ --> kg
    T_0 = 34 + 273; %@@@ --> K
    T_1 = 89 + 273; %@@@ --> K
    c = 4186; % J/(kg K)

    dS = m * c * log(T_1/T_0);
    fprintf("dS=  %0.3f J/K\n", dS);
end


%%%%%~~~~
if select == 14
    V = 0.019; %@@@ m^3
    T_0 = 22 + 273; %@@@ --> K
    P_0 = 9e5; % Pa
    
    df = 3 + 2 + 2;
    fprintf("a) df=   %d\n", df);
    fprintf("b) isochroric ... volume don't change\n");
    dT = 10; % temperature increases
    T_1 = T_0 + dT;
    P = (T_1) * P_0 / T_0;
    fprintf("c) P =  %0.1f  Pa\n", P);

    dU = (df/2) * (P_0 * V / T_0) * dT;
    fprintf("d) dU=  %0.3f J\n", dU);

    dS = (df/2) * (P_0 * V / T_0) * log(T_1/T_0);
    fprintf("e) dS=  %0.3f J/K\n", dS);

end


%%%%%~~~~
if select == 99
    fprintf("\n\tDONE\n");
end


%%%%%%%%~~~~~~~~END>  v2_ch3_ch4_thermod.m
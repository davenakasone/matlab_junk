%{
    ee442, hw3

    1 : p2
    2 : p3
    3 : p4
%}
clc;
clear;
select = 3;


if select == 1
    R        = 1e3;        % ohms, from resistor
    V_F      = 10;         % V, from pulse source
    V_R      = -10;        % V, from pulse source
    VD_sil   = 0.55;       % V, from model
    VD_shot  = 0.2;        % V, from model
    Cj0_sil  = 45e-12;     % F, from model
    Cj0_shot = 190e-12;    % F, from model
    tau_sil  = 45e-9;      % s, from model
    tau_shot = 0;          % s, from model

    t_1     = 10e-9;        % time when the pulse switches
    t2_sil  = 37.85e-9;     % when silicon diode is done storage
    t2_shot = t_1;          % when the schottky diode is done storage
    t3_sil  = 91.06e-9;     % when the silicon diode is recovered
    t3_shot = 354.14e-9;    % when the schottcky diode is recovered

    % find forward and reverse currents
    IF_sil  = (V_F - VD_sil) / R;
    IR_sil  = (V_R - VD_sil) / R;
    IF_shot = (V_F - VD_shot) / R;
    IR_shot = (V_R - VD_shot) / R;

    % compare calculated and observed storage times
    ts_sil  = tau_sil  * log((IF_sil - IR_sil)   / (-1 * IR_sil));
    ts_shot = tau_shot * log((IF_shot - IR_shot) / (-1 * IR_shot));
    fprintf("theoretical, silicon diode , t_store =  %0.3f ns\n", ts_sil * 1e9);
    fprintf("actual     , silicon diode , t_store =  %0.3f ns\n", (t2_sil-t_1) * 1e9);
    fprintf("theoretical, schottky diode, t_store =  %0.3f ns\n", ts_shot * 1e9);
    fprintf("actual     , schottky diode, t_store =  %0.3f ns\n", (t2_shot-t_1) * 1e9);

     % compare calcuated and observed recovery times
     tr_sil  = t3_sil - t_1;
     tr_shot = t3_shot - t_1;
     fprintf("\nsilicon diode , t_rr = %0.3f ns\n", tr_sil * 1e9);
     fprintf("schottky diode, t_rr = %0.3f ns\n", tr_shot * 1e9);

     % determine the capacitance, VD = 10% of original
     C_sil  = tr_sil / (-1 * R * log(0.1));
     C_shot = tr_shot / (-1 * R * log(0.1));
     fprintf("\nsilicon Cj0  =  %7.3f pF  ,  Cj = %7.3f pF\n", Cj0_sil*1e12, C_sil*1e12);
     fprintf("schottky Cj0 =  %7.3f pF  ,  Cj = %7.3f pF\n", Cj0_shot*1e12, C_shot*1e12);

end


%%%%~~~~~


if select == 2
    V_m = 120 * sqrt(2);
    R = 18;

    V_o_avg = 2 * V_m / pi;
    I_o_avg = V_o_avg / R;
    fprintf("average output current, I_o_avg  =  %0.3f A\n", I_o_avg);
    I_o_peak = V_m / R;
    fprintf("peak output current   , I_o_peak =  %0.3f A\n", I_o_peak);
    I_o_rms = I_o_peak / sqrt(2);
    fprintf("rms output current    , I_o_rms  =  %0.3f A\n", I_o_rms);
    I_d_avg = I_o_avg / 2;
    fprintf("\naverage diode current , I_d_avg  =  %0.3f A\n", I_d_avg);
    I_d_peak = I_o_peak;
    fprintf("peak diode current    , I_d_peak =  %0.3f A    == I_o_peak\n", I_d_peak);
    I_d_rms = I_o_rms / sqrt(2);
    fprintf("rms diode current     , I_d_rms  =  %0.3f A\n", I_d_rms);
end


%%%%~~~~~


if select == 3
    f = 60;
    w = 2 * pi * f;
    V_d = 0.55;
    V_in_rms = 120;
    V_in_peak = V_in_rms * sqrt(2);
    L_p = 1e3;
    L_s = 2.75;
    V_out_peak = V_in_peak / sqrt(L_p/L_s)
    V_out_rms = V_out_peak / sqrt(2);
    I_out_max = 100e-3;
    R_max = V_out_peak / I_out_max;
    V_rip_max = 10e-3;
    V_rip_target = V_rip_max / 10;

    
    dots = 10;

%     theta = atan(-1 * w * R * C);
%     syms alp;
%     eqn_a = sin(alp);
%     eqn_b = sin(theta) * exp( (-1 * (2*pi + alp - theta) ) / (w * R * C));
%     alpha = vpasolve(eqn_a == eqn_b, alp);
%     check = double(subs(eqn_a - eqn_b, alp, alpha));
%     fprintf("\nchecked:  %0.3f, should be 0\n", check);
%     vripple_exact = V_m * (1 - sin(-1*alpha));
%     fprintf("b) exact peak-peak ripple voltage:  %0.4f V\n", vripple_exact);
end

%%%%~~~~~


if select == 99
end


%%%%%%%%~~~~~~~~END>  hw3.m

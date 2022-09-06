%{
    chapter 7 : steady state error

    1  :  ex7.1, p624 basic error calc s->0
    2  :  ex7.2, p627 steady state error, no integrators
    3  :  ex7.3, p628 steady state error, 1 integrator in path
    4  :  sa7.1, p629 basic steady state error
    5  :  ex7.4, p631 static error const, Kp, Kv, Ka
    6  :  sa7.2, unity fedback
    7  :  ex7.6, only ramp has finite error in type 1
    8  :  sa7.3, type0, find K
    9  :  hw7.1, find error
    10 :  hw7.10, find unknown in a limit
    11 :  hw7.11, find K, see the parabola
    12 :  hw7.12, reduced to unity feedback
    13 :  hw7.13, reduce to unity feedback
    14 :  hw7.19, system type, K unknown
    15 :  hw7.25
    16 :  hw7.26

%}
close all;
clc;
select = 16;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    Ts = 5 / (s^2 + 7*s +10);
    Rs = 1 / s; % unit step
    Es = Rs * (1 - Ts);
    ss_error = limit(s * Es, 0)
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    Gs = 120 * (s + 2) / ((s + 3) * (s + 4));
    Rs_a = 5/s;
    Rs_b = 5/s^2;
    Rs_c = 5/s^3;
    
    e_a = limit((s * Rs_a) / (1 + Gs), s, 0)
    e_b = limit((s * Rs_b) / (1 + Gs), s, 0)
    e_c = limit((s * Rs_c) / (1 + Gs), s, 0)
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    Gs = 100 * (s + 2) * (s + 6) / (s * (s + 3) * (s + 4));
    Rs_a = 5/s;
    Rs_b = 5/s^2;
    Rs_c = 5/s^3;
    
    e_a = limit((s * Rs_a) / (1 + Gs), s, 0)
    e_b = limit((s * Rs_b) / (1 + Gs), s, 0)
    e_c = limit((s * Rs_c) / (1 + Gs), s, 0)
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    Gs = 10 * (s + 20) * (s + 30) / (s * (s + 25) * (s + 35));
    Rs_a = 15/s;
    Rs_b = 15/s^2;
    Rs_c = 15/s^3;
    e_a = limit((s * Rs_a) / (1 + Gs), s, 0)
    e_b = limit((s * Rs_b) / (1 + Gs), s, 0)
    e_c = limit((s * Rs_c) / (1 + Gs), s, 0)
    Gs = 10 * (s + 20) * (s + 30) / (s^2 * (s + 25) * (s + 35) * (s+50));
    pretty(simplify(expand(Gs)))
    rootz = roots([1, 110, 3875, 43750, 0, 0])  % unstable, can't find
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    R_step = 1 / s;
    R_ramp = 1 / s^2;
    R_para = 1 / s^3;
    Gs_a = (500 * (s + 2) * (s + 5)) / ((s + 8) * (s + 10) * (s + 12));
    Gs_b = (500 * (s + 2) * (s + 5) * (s + 6)) / (s * (s + 8) * (s + 10) * (s + 12));
    Gs_c = (500 * (s + 2) * (s + 4) * (s + 5) * (s + 6) * (s + 7)) / (s^2 * (s + 8) * (s + 10) * (s + 12));
    
    Kp_a = double(limit(Gs_a, s, 0));
    Kv_a = double(limit(s * Gs_a, s, 0));
    Ka_a = double(limit(s^2 * Gs_a, s, 0));
    fprintf("Gs_a>>  Kp=  %0.3f  ,  Kv=  %0.3f  ,  Ka=  %0.3f\n", Kp_a, Kv_a, Ka_a);
    Kp_b = double(limit(Gs_b, s, 0));
    Kv_b = double(limit(s * Gs_b, s, 0));
    Ka_b = double(limit(s^2 * Gs_b, s, 0));
    fprintf("Gs_b>>  Kp=  %0.3f  ,  Kv=  %0.3f  ,  Ka=  %0.3f\n", Kp_b, Kv_b, Ka_b);
    Kp_c = double(limit(Gs_c, s, 0));
    Kv_c = double(limit(s * Gs_c, s, 0));
    Ka_c = double(limit(s^2 * Gs_c, s, 0));
    fprintf("Gs_c>>  Kp=  %0.3f  ,  Kv=  %0.3f  ,  Ka=  %0.3f\n", Kp_c, Kv_c, Ka_c);
    
    e_ss_step_a = 1 / (1 + Kp_a);
    e_ss_ramp_a = 1 / Kv_a;
    e_ss_para_a = 1 / Ka_a;
    fprintf("\nGs_a>> e_step:  %0.3f  ,  e_ramp:  %0.3f  ,  e_para:  %0.3f\n",...
        e_ss_step_a, e_ss_ramp_a, e_ss_para_a);
    e_ss_step_b = 1 / (1 + Kp_b);
    e_ss_ramp_b = 1 / Kv_b;
    e_ss_para_b = 1 / Ka_b;
    fprintf("Gs_b>> e_step:  %0.3f  ,  e_ramp:  %0.3f  ,  e_para:  %0.3f\n",...
        e_ss_step_b, e_ss_ramp_b, e_ss_para_b);
    e_ss_step_c = 1 / (1 + Kp_c);
    e_ss_ramp_c = 1 / Kv_c;
    e_ss_para_c = 1 / Ka_c;
    fprintf("Gs_b>> e_step:  %0.3f  ,  e_ramp:  %0.3f  ,  e_para:  %0.3f\n",...
        e_ss_step_c, e_ss_ramp_c, e_ss_para_c);
end


%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    Gs = (1000 * (s + 8)) / ((s + 7) * (s +9));
    Ts = Gs / (1 + Gs);
    
    Kp = double(limit(Gs, s, 0))
    Kv = double(limit(s * Gs, s, 0))
    Ka = double(limit(s^2 * Gs, s, 0))
    
    e_step = 1 / (1 + Kp)
    e_ramp = 1 / Kv
    e_para = 1 / Ka
end


%------------------------------------------------------------------------------------------
if (select == 7)
    % wants error = 10%    e = 1/Kv = 0.1, Kv = 10
    syms K;
    eqn = (K * 5 ) / (6 * 7 * 8);
    k = solve(eqn == 10, K) 
end


%------------------------------------------------------------------------------------------
if (select == 8)
    syms K;
    eqn = 1 / (1 +((K * 12) / (14 * 18)));
    k = solve(eqn==0.1, K);
end

%------------------------------------------------------------------------------------------
if (select == 9)
    syms s;
    Gs = (1350 * (s + 2) * (s + 10) * (s + 32)) / (s * (s + 4) * (s^2 + 8*s + 32));
    Rs_a = 17/s;
    Rs_b = 32/s^2;
    Rs_c = 48/s^3;
    
    ss_err_a = limit((s * Rs_a) / (1 + Gs), s, 0);
    ss_err_b = limit((s * Rs_b) / (1 + Gs), s, 0);
    ss_err_c = limit((s * Rs_c) / (1 + Gs), s, 0);
    fprintf("\nsteady state error, R(s) = 17 u(t)::  %0.6f\n", ss_err_a);
    fprintf("\nsteady state error, R(s) = 32 t u(t)::  %0.6f\n", ss_err_b);
    fprintf("\nsteady state error, R(s) = 48 t^2 u(t)::  %0.6f\n", ss_err_c);
    
end


%------------------------------------------------------------------------------------------
if (select == 10)
    syms alpha;
    syms s;
    Kv = 40e3;
    Gs = (300e3 * (s + 5) * (s + 10) * (s + 30)) /...
        (s * (s + 60) * (s + alpha) * (s + 90));
    eqn = limit(s * Gs, s, 0);
    aa = solve(eqn == Kv, alpha);
    fprintf("\nalpha = %0.3f\n", aa);
end


%------------------------------------------------------------------------------------------
if (select == 11)
    syms s;
    syms K;
    K_a = 10e3;
    Gs = (K * (s + 2) * (s + 4) * (s + 6)) /...
        (s^2 * (s + 5) * (s + 7));
    eqn = limit( s^2 * Gs, s, 0);
    k = solve(eqn == K_a, K);
    fprintf("K should be:  %0.3f\n", k);
end


%------------------------------------------------------------------------------------------
if (select == 12)
    syms s;
    Rs_step = 20 / s;
    Rs_ramp = 20 / s^2;
    Rs_para = 20 / s^3;
    Gs = 6 / (s * (s +1) * (s + 5));
    Hs = s + 2;
    Gs_p = Gs / (1 + Gs * Hs);
    %pretty(expand(simplify(Gs_p)));
    
    Kp = limit(Gs_p, s, 0);
    Kv = limit(s * Gs_p, s, 0);
    Ka = limit(s^2 * Gs_p, s, 0);
    fprintf("Kp =  %0.3f\n", Kp);
    fprintf("Kv =  %0.3f\n", Kv);
    fprintf("Ka =  %0.3f\n", Ka);
    
    ss_err_step = limit(s * Rs_step / (1 + Gs_p), s, 0);
    ss_err_ramp = limit(s * Rs_ramp / (1 + Gs_p), s, 0);
    ss_err_para = limit(s * Rs_para / (1 + Gs_p), s, 0);
    fprintf("steady state error, step:  %0.3f\n", ss_err_step);
    fprintf("steady state error, ramp:  %0.3f\n", ss_err_ramp);
    fprintf("steady state error, para:  %0.3f\n", ss_err_para);
    
end


%------------------------------------------------------------------------------------------
if (select == 13)
    syms s;
    Gs_a = (100 * (s + 2)) / (s * (s + 5));
    Gs_b = 1000 / s;
    Ms = 5;
    Ts = (Gs_a / (1 + Gs_a * Ms)) * Gs_b;
    pretty(simplify(expand(Ts)));
end


%------------------------------------------------------------------------------------------
if (select == 14)
    syms s;
    syms K;
    Gs = (K * (s^2 + 6*s + 6)) / ((s + 5)^2 * (s + 3));
    %pretty(expand(Gs));
    
    Rs_a = 12 / s;
    Rs_b = 12 / s^2;
    ss_err_a = limit(s * Rs_a / (1 + Gs), s, 0);
    ss_err_b = limit(s * Rs_b / (1 + Gs), s, 0);
    pretty(simplify(expand(ss_err_a)));
    pretty(simplify(subs(ss_err_b, K,1)));
end


%------------------------------------------------------------------------------------------
if (select == 15)
    o_s = 0.1;
    zeta = -log(o_s) / sqrt(pi^2 + (log(o_s)^2))
end


%------------------------------------------------------------------------------------------
if (select == 16)
    syms s;
    Gs_a = 1 / (s^2 * (s + 1));
    Gs_b = 1 / (s^2 * (s + 3));
    Hs = 1 / s;
    Ts_p = Gs_a / (1 + Gs_a * Hs);
    Gs = Ts_p * Gs_b;
    %fprintf("G(s)\n");
    %pretty(simplify(Gs));
    Ts = Gs / (1 + Gs);
    %fprintf("T(s)\n");
    %pretty(factor(Ts));
    %rotz = double(solve(1/Ts==0, s))
    err_s = s / (1 + Gs);
    pretty(simplify(err_s))
    inp_a = 5/s;
    inp_b = 5/s^2;
    err_a = limit(inp_a * err_s, s, 0)
    err_b = limit(inp_b * err_s, s, 0)
    
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end
%%%%%%%%~~~~~~~~END>
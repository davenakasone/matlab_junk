%{
    chapter 11 : frequency design

    1  :  ex11.1, design
    2  :  sa11.1, design
    3  :  ex11.2, the lag compensator design

    rlocus()
    rlocusplot()
    sgrid()
    step()
    Control System Designer
    nyquist()  ...right click and say all margins
    nicholas()
    
    margin()
    APPS --> control --> linear system analyzer
    controlSystemDesigner
%}
format compact;
clear;
close all;
clc;
select = 3;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    syms w;
    o_s = 0.095;
    K = 3.6;
    Ls = K * 100 / (s * (s + 100) * (s + 36));
    
    tf_Ls = f_Fs_2_tf(Ls);
    zeta = f_os_2_zeta(o_s);
    p_m = f_zeta_2_pm_rad(zeta);
    fprintf("this implies zeta=  %0.3f  ,  Phi_m=  %0.3f deg\n", zeta, rad2deg(p_m));
    p_target = -1 * (pi - p_m);
    fprintf("taking difference from -180, target angle is:  %0.3f  deg\n", rad2deg(p_target));
    w_pm = 14.8;
    fprintf("bode plot says w = %0.3f rad/sec where phase margin occurs\n", w_pm);
    M_pm = subs(abs(Ls), s, 1j*w_pm);
    M_pm_db = 20*log10(M_pm);
    fprintf("at w_pm, M=  %0.3f\n", M_pm_db);
    increase = 10^(-1*M_pm_db/20)
    fprintf("must increase:  %0.3f  , K = %0.1f\n", double(increase), increase*K);
    bode(tf_Ls, "g-");
    nyquist(tf_Ls);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    syms w;
    o_s = 0.2;
    K = 1;%194200;
    Ls = K / (s * (s + 50) * (s + 120));
    
    tf_Ls = f_Fs_2_tf(Ls);
    zeta = f_os_2_zeta(o_s);
    p_m = f_zeta_2_pm_rad(zeta);
    fprintf("this implies zeta=  %0.3f  ,  Phi_m=  %0.3f deg\n", zeta, rad2deg(p_m));
    p_target = -1 * (pi - p_m);
    fprintf("taking difference from -180, target angle is:  %0.3f  deg\n", rad2deg(p_target));
    w_pm = 27.4;
    fprintf("bode plot says w = %0.3f rad/sec where phase margin occurs\n", w_pm);
    M_pm = subs(abs(Ls), s, 1j*w_pm);
    M_pm_db = 20*log10(M_pm);
    fprintf("at w_pm, M=  %0.3f\n", M_pm_db);
    increase = 10^(-1*M_pm_db/20)
    fprintf("must increase:  %0.3f  ,  K=  %0.1f\n", double(increase), increase*K);
    %bode(tf_Ls, "g-");
    %nyquist(tf_Ls);
    controlSystemDesigner(tf_Ls)
end


%------------------------------------------------------------------------------------------
if (select == 3)
    % G_lag = (s + 1/T) / (s + 1/aT)  , a > 1
    % should attenuate low freq
    % 4 steps, p1112
    syms s;
    syms w;
    juice = deg2rad(12);
    o_s = 0.095;
    improve = 10; % improve Kv this much
    K = 5839; % you want a 10x error improvement from the lead compensated design in ex11.1
    Ls = K / (s * (s + 36) * (s + 100));

    zeta = f_os_2_zeta(o_s);
    p_m = f_zeta_2_pm_rad(zeta);
    p_target = -1 * (pi - p_m + juice);
    Kv = limit(s*Ls, s, 0);
    Kv_new = Kv*improve;

    fprintf("give OS=  %0.3f  ,  zeta=  %0.3f\n", o_s, zeta);
    fprintf("this means phase margin=  %0.3f deg\n", rad2deg(p_m));
    fprintf("makes target=  %0.3f deg\n", p_target);
    fprintf("Kv=  %0.3f  ,  improve to Kv=  %0.3f\n", Kv, Kv_new);
    
    
    

end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
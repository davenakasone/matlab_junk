%{
    quiz3, ch10, ch11

    1  :  p1
    2  :  p2
    3  :  check David Pinnale's solution for p2


    rlocus()
    rlocusplot()
    sgrid
    ngrid
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
select = 2;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    syms w;
    K = 100;
    Ls = K / (s * (s + 5)^2);
    tf_Ls = f_Fs_2_tf(Ls);

    %{
    figure;
    hold on;
    axis equal;
    nyquist(tf_Ls, "b-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;
    %}
    %
    figure;
    hold on;
    bode(tf_Ls, "r-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    syms w;
    K = 2000;
    Ls = K / (s * (s + 10)^2);
    tf_Ls = f_Fs_2_tf(Ls);
    zeta = 0.6;
    p_m = f_zeta_2_pm_rad(zeta);
    fprintf("phase margin=  %0.3f deg, or approx=  %0.3f deg\n", rad2deg(p_m), 100*zeta);
    p_target = -1 * (pi - (p_m + deg2rad(10)));
    fprintf("to be safe, find when the phase is:  %0.3f\n", rad2deg(p_target));
    %{
    figure;
    hold on;
    bode(tf_Ls, "c-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;
    %}
    w_c = 1.75;
    db_ucomp = 20.9;
    fprintf("\nat %0.3f deg,  w_c=  %0.3f  rad/sec ,  M=  %0.3f dB\n",...
        rad2deg(p_target), w_c, db_ucomp);
    Ljw = subs(Ls, s, 1j*w_c);
    check_wc = angle(Ljw);
    check_M = 20*log10(abs(Ljw));
    fprintf("angle check=  %0.3f deg  ,  M check=  %0.3f dB\n", rad2deg(check_wc), check_M);
    
    z_lag = w_c / 10;
    tau = 1 / z_lag;
    alpha = 10^(db_ucomp/20);
    fprintf("\ntherefore, the zero is one decade below, w_c/10 =  %0.3f\n", z_lag);
    fprintf("solving for alpha = %0.3f\n", alpha);
    p_lag = 1 / (alpha * tau);
    fprintf("compensator zero:  %0.3f  [%0.3f],  compensator pole:  %0.3f [%0.3f]\n",...
        -1*z_lag, 1/z_lag, -1*p_lag, 1/p_lag);
    G_lag = (1 + s*(1/z_lag)) / (1 + s*(1/p_lag));
    tf_G_lag = f_Fs_2_tf(G_lag);
    G_comp = Ls * G_lag;
    tf_comp = f_Fs_2_tf(G_comp);
    Kv_now = limit(s*G_comp, s, 0);

    %{
    figure;
    hold on;
    axis equal;
    nyquist(tf_Ls, "r--");
    nyquist(tf_G_lag, "b--");
    nyquist(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
    %{
    figure;
    hold on;
    bode(tf_Ls, "r--");
    bode(tf_G_lag, "b--");
    bode(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    syms w;
    K = 2000;
    Ls = K / (s * (s + 10)^2);
    tf_Ls = f_Fs_2_tf(Ls);
    zeta = 0.6;
    p_m = f_zeta_2_pm_rad(zeta);
    fprintf("phase margin=  %0.3f deg\n", rad2deg(p_m));

    var_n = 4.367;
    var_d = 35.71;
    G_lag = (1 + s*var_n) / (1 + s*var_d );
    tf_G_lag = f_Fs_2_tf(G_lag);
    G_comp = Ls * G_lag;
    tf_comp = f_Fs_2_tf(G_comp);

    %
    figure;
    hold on;
    axis equal;
    nyquist(tf_Ls, "r--");
    nyquist(tf_G_lag, "b--");
    nyquist(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
    %
    figure;
    hold on;
    bode(tf_Ls, "r--");
    bode(tf_G_lag, "b--");
    bode(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
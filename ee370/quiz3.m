%{
    quiz3, ch10, ch11

    1  :  p1
    2  : p2


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
    fprintf("phase margin=  %0.3f deg\n", rad2deg(p_m));

    tau = 0.2;
    alpha = 2;
    G_lag = (s + (1/tau)) / (s + (1/(alpha*tau)));
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
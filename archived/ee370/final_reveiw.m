%{
    final review

    1  :  Routh, special case 1, 0 in a column
    2  :  Routh, special case 2, row becomes 0
    3  :  frequency
    4  : Gm, not on test
    5  : get phase margin
    6  : design

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
select = 6;




%------------------------------------------------------------------------------------------
if (select == 1)
    % given char: s^5 + 2*s^4 + 3*s^3 + 6*s^2 + 5*s + 3
    syms s;
    syms ep;
    a_3_1 = f_routh(1,3,2,6,2);
    a_3_1 = ep; % because it was 0
    a_3_2 = f_routh(1, 5, 2, 3, 2);
    
    Ts = 1/ (s^5 + 2*s^4 + 3*s^3 + 6*s^2 + 5*s + 3);
    tf_Ts = f_Fs_2_tf(Ts);
    pole(tf_Ts)
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    Ts = 1 / (s^5 + 7*s^4 + 6*s^3 + 42*s^2 + 8*s + 56);
    tf_Ts = f_Fs_2_tf(Ts);
    pole(tf_Ts)
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    syms w;
    Gs = 100 * (s + 5) / ((s + 2) * (s + 9) * s);
    tf_Gs = f_Fs_2_tf(Gs);

    
    
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    syms w;
    Gs = 1/ (s * (s + 2) * (s + 3));
    tf_Gs = f_Fs_2_tf(Gs);
    Gw_angle = 0 - ((pi/2) + atan(w/2) + atan(w/3));
    Gw_abs = 1/ (w * sqrt(w^2 + 4) * sqrt(w^2 + 9));

    nyquist(tf_Gs);
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    syms w;
    K = 100;
    Gs = K / (s* (s + 5));
    tf_Gs = f_Fs_2_tf(Gs);

    %nyquist(tf_Gs);
    fprintf("");
end


%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    syms w;
    Ls = 10 / (s * (s + 2));
    Lw_abs = 10 / (w * sqrt(w^2 + 4));
    tf_Ls = f_Fs_2_tf(Ls);
    zeta = 0.6;
    p_m = f_zeta_2_pm_rad(zeta);
    p_target = -1*pi + p_m + deg2rad(5);
    fprintf("phase margin should be:  %0.3f deg\n", rad2deg(p_m));

    
    %{
    figure;
    hold on;
    axis equal;
    nyquist(tf_Ls, "r--");
    %nyquist(tf_G_lag, "b--");
    %nyquist(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
    %{
    figure;
    hold on;
    bode(tf_Ls, "r--");
    %bode(tf_G_lag, "b--");
    %bode(tf_comp, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("uncomp", "lag", "COMP", location='best', fontsize=14);
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
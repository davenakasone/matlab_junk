%{
    chapter 9, homework

    p4  :  lag design using ex9.2
    p12 :  lead design using ex9.4
    p15 :  lag-lead design using ex9.6
    p20 :  PID design, using ex9.5

%}
format compact;
clear;
close all;
clc;
select = 4;


%------------------------------------------------------------------------------------------
if (select == 1)
    os_initial = 0.1;
    zeta_initial = f_os_2_zeta(os_initial);
    fprintf("\ngiven the OS:  %0.1f %%  ,  zeta = %0.3f\n",...
        100*os_initial, zeta_initial);

    syms s;
    Gs = 1 / ((s + 2) * (s + 4) * (s + 6));
    tf_initial = f_Fs_2_tf(Gs);
    %f_rlocus(tf_initial);

    K = 45.4;
    p_new = -2.03 + 1j*2.76;
    K_p = limit(K*Gs, s, 0);
    Kp_improve = 20;
    fprintf("\nKp = %0.3f  ,  uncompensated\n", K_p);
    fprintf("to improve to Kp = %d  ,  must increase by:  %0.3f\n",...
        Kp_improve, Kp_improve/K_p);

    p_lag = 0.01;
    z_lag = (Kp_improve / K_p) * p_lag;
    fprintf("use G_lag(s) = (s + %0.3f) / (s + %0.3f)\n", z_lag, p_lag);
    Gs_lag = (s + z_lag) / (s + p_lag);
    Gs_improved = Gs * Gs_lag;
    tf_improved = f_Fs_2_tf(Gs_improved);
    %f_rlocus(tf_improved);

    K_improved = 45.8;
    Kp_new = limit(K_improved * Gs_lag * Gs, s, 0);
    fprintf("\nKp = %0.3f  ,  when compensated, the goal is achieved\n", Kp_new);

    tu = K * Gs / (1 + K * Gs);
    tf_ucomp = f_Fs_2_tf(tu);
    tc = K_improved * Gs_lag * Gs / (1 + K_improved * Gs_lag * Gs);
    tf_comp = f_Fs_2_tf(tc);

    t_final = 100;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(tf_ucomp, t_final, "r-");
    step(tf_comp, t_final, "b-");
    grid on;
    legend({'Uncompensated', 'Compensated'},...
        Location='best', FontSize=14);
    templ = findall(gcf,'type','line');
    set([templ(5), templ(7)], 'linewidth',2);
    set([templ(1), templ(2)], 'linewidth',1);
    %set(groot,'defaultLineMarkerSize',15);
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 2)
    os_given = 0.205;
    zeta = f_os_2_zeta(os_given);
    fprintf("\ngiven the OS:  %0.1f %%  ,  zeta = %0.3f\n",...
        100*os_given, zeta);
    Ts_given  = 3;
    w_n = 4 / (zeta * Ts_given);
    fprintf("\ngiven Ts = %d  ,  w_n = %0.3f\n", Ts_given, w_n);
    
    syms s;
    Gs = 1 / (s^2 * (s + 4) * (s + 12));
    tf_initial = f_Fs_2_tf(Gs);
    %f_rlocus(tf_initial);

    s_dom1 = -1 * zeta * w_n + 1j * w_n * sqrt(1 - zeta^2);
    s_dom2 = conj(s_dom1);
    fprintf("\nthe compensated dominate poles:\n");
    fprintf("s = %0.3f + j%0.3f\n", real(s_dom1), imag(s_dom1));
    fprintf("s = %0.3f - j%0.3f\n", real(s_dom2), -1*imag(s_dom2));
    
    rootz = roots([1, 16, 48, 0, 0])
    z_lead = .01;
    vzlead = angle(s_dom1 - (-1)*z_lead);
    v1 = angle(s_dom1 - rootz(1));
    v2 = angle(s_dom1 - rootz(2));
    v3 = angle(s_dom1 - rootz(3));
    v4 = angle(s_dom1 - rootz(4));
    angle_total = vzlead - (v1 + v2 + v3 + v4);
    fprintf("angles without pole sum to:  %0.3f degree\n", rad2deg(angle_total));
    plead_angle = pi + angle_total;
    fprintf("angle of pole must be:  %0.3f degree\n", rad2deg(plead_angle));
    p_lead = -1*real(s_dom1) + (imag(s_dom1)/tan(plead_angle));
    fprintf("p_lead = %0.3f\n", p_lead);

    Gs_lead = (s + z_lead) / (s + p_lead);
    tf_improved = f_Fs_2_tf(Gs * Gs_lead);
    %f_rlocus(tf_improved);

    K = 4.21e3;
    G_imp = K * Gs_lead * Gs;
    T_imp = G_imp / (1 + G_imp);
    tf_imp = f_Fs_2_tf(T_imp);
    tf_zeros = zero(tf_imp)
    tf_poles = pole(tf_imp)
    figure();
    pzplot(tf_imp);
    ttt = findobj(gca, 'type', 'line');
    set(ttt, 'markersize', 10);
    set(ttt, 'linewidth', 2);
    hold on;
    grid on;
    axis padded;
    hold off;

    t_final = 10;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(tf_imp, t_final, "r-");
    grid on;
    templ = findall(gcf,'type','line')
    set([templ(4), templ(3)], 'linewidth',3);
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    temp_u = (355)/((s+2)*(s+4)*(s+6)*(s+8));
    u_tf = f_Fs_2_tf(temp_u/(1+temp_u));
    temp_l = (1360*(s+5))/((s+2)*(s+4)*(s+6)*(s+8)*(s+14.679));
    l_tf = f_Fs_2_tf(temp_l/(1+temp_l));
    temp_c = (1417*(s+5)*(s+0.04759))/((s+2)*(s+4)*(s+6)*(s+8)*(s+0.001)*(s+14.679));
    c_tf = f_Fs_2_tf(temp_c/(1+temp_c));
    t_final = 150;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(u_tf, t_final, "b-");
    step(l_tf, t_final, "g-");
    step(c_tf, t_final, "r-");
    grid on;
    legend({'Uncompensated', 'Lead Only' ,'Compensated'},...
        Location='best', FontSize=14);
    temp = findall(gcf,'type','line')
    set([temp(6), temp(7), temp(8), temp(9), temp(10)], 'linewidth',2);
    set([temp(1), temp(2)], 'linewidth',1);
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    Gs = 1 /((s + 4) * (s + 6) * (s + 10));

    o_s = 0.25;
    Ts_limit = 2;
    zeta_wn = 4/Ts_limit;
    fprintf("\nif limited to Ts= 2, zeta * w_n = %0.1f\n", zeta_wn);
    zeta = f_os_2_zeta(o_s);
    theta = acos(zeta);
    fprintf("from the %%OS= %0.1f  ,  zeta=  %0.3f  ,  theta=  %0.3f degree\n",...
        100*o_s, zeta, rad2deg(theta));
    w_n = 2 / zeta;
    fprintf("making w_n =  %0.3f\n", w_n);
    op_r = -1 * zeta_wn;
    op_i = 1j* w_n * sqrt(1-zeta^2);
    s_op = op_r + op_i;
    fprintf("this makes the operating point:  %0.3f + j%0.3f\n", real(s_op), imag(s_op));

    p_c = 0;
    p_1 = -4;
    p_2 = -6;
    p_3 = -10;
    v_a = angle(s_op - p_c);
    v_b = angle(s_op - p_1);
    v_c = angle(s_op - p_2);
    v_d = angle(s_op - p_3);
    v_total = 2*pi-(v_a + v_b + v_c + v_d);
    fprintf("\nsum of angles without zero:  %0.3f degree\n", rad2deg(v_total));
    theta_zc = pi-v_total;
    fprintf("the angle of the zero on the PD should be:  %0.3f\n", rad2deg(theta_zc));
    z_c = (imag(s_op) / tan(theta_zc)) - real(s_op);
    fprintf("the zero of the PD s=  %0.3f\n", -1*z_c);

    Gs_pd = (s + z_c) / s;
    tf_pd = f_Fs_2_tf(Gs * Gs_pd);
    %f_rlocus(tf_pd);
    K = 294;
    Gs_pi = (s + 0.001) / s;
    Gs_pid = K * Gs * Gs_pd * Gs_pi;
    Ts_pid = Gs_pid / (1 + Gs_pid);
    tf_pid = f_Fs_2_tf(Ts_pid);

    t_final = 5;
    figure(Position=[20, 20, 800, 800]);
    step(tf_pid, t_final, "b-");
    hold on;
    grid on;
    title("step response of the PID", FontSize=20);
    temp = findall(gcf,'type','line')
    set([temp(3), temp(4)], 'linewidth',2);
    set([temp(1), temp(2)], 'linewidth',1);
    hold off;

    Hs_ramp = (1/s) * Ts_pid;
    tf_ramp = f_Fs_2_tf(Hs_ramp);
    t_final = 50;
    figure(Position=[20, 20, 800, 800]);
    step(tf_ramp, t_final, "r-");
    hold on;
    grid on;
    title("ramp response of the PID", FontSize=20);
    temp = findall(gcf,'type','line')
    set([temp(3), temp(4)], 'linewidth',2);
    set([temp(1), temp(2)], 'linewidth',1);
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>

%{
    zeta = 0.5;
    theta = acos(zeta);
    os_given = f_zeta_2_os(zeta);
    fprintf("\ngiven zeta, theta=  %0.1f degree,  the %%OS=  %0.3f\n",...
        rad2deg(theta), 100*os_given);

    syms s;
    Gs = 1 / ((s + 2) * (s + 4) * (s + 6) * (s + 8));
    tf_initial = f_Fs_2_tf(Gs);
    %f_rlocus(tf_initial);

    K_u = 355;
    Gs_u = K_u * Gs;
    tf_u = f_Fs_2_tf(Gs_u);
    Kp_u = limit(Gs_u, s, 0);
    fprintf("Kp = %0.3f\n", Kp_u);
    e_u = 1 / (1 + Kp_u);
    fprintf("uncompensated error:  %0.3f\n", e_u);

    s_operate = -1.53 + 1j * 2.65;
    Ts = 4 / abs(real(s_operate));
    fprintf("\nthe current Ts = %0.3f seconds\n", Ts);
    Ts_new = Ts - 0.5;
    fprintf("the design needs Ts = %0.3f\n", Ts_new);
    zeta_wn = -1 * 4 / Ts_new;
    fprintf("this means zeta * wn =  %0.3f\n", zeta_wn);
    s_op_i = -1 * tan(theta) * zeta_wn;
    s_op = zeta_wn + 1j * s_op_i;
    fprintf("operate now at s= %0.3f + j%0.3f\n",...
        real(s_op), imag(s_op));
    
    z_lead = 5;
    v_z = angle(s_op - (-1)*z_lead);
    v_1 = angle(s_op - (-1)*2);
    v_2 = angle(s_op - (-1)*4);
    v_3 = angle(s_op - (-1)*6);
    v_4 = angle(s_op - (-1)*8);
    anglez = v_z - (v_1 + v_2 + v_3 + v_4);
    fprintf("without the compensator pole, the angles are:  %0.3f degree\n", rad2deg(anglez));
    pc_lead_ang = pi + anglez;
    fprintf("the vector angle of the pole:  %0.3f degree\n", rad2deg(pc_lead_ang));
    p_lead = -1*real(s_op) + (imag(s_op)/tan(pc_lead_ang));
    fprintf("the poles is at:  %0.3f\n", -1*p_lead);

    Gs_lead = (s + z_lead) / (s + p_lead);
    Gss = Gs * Gs_lead;
    tf_lead = f_Fs_2_tf(Gss);
    %f_rlocus(tf_lead);

    K_l = 1360;
    Gss_lead = K_l * Gss;
    Kp_l = limit(Gss_lead, s, 0);
    fprintf("Kp = %0.3f\n", Kp_l);
    e_l = 1 / (1 + Kp_l);
    fprintf("with lead error:  %0.3f\n", e_l);
    e_target = e_u/30;
    fprintf("make the error:  %0.3f\n", e_target);
    improve = e_l / e_target;
    fprintf("improve this much:  %0.3f\n", improve);
    e_ll = e_l / improve;
    fprintf("final error is:  %0.3f\n", e_ll);
    Kp_ll = (improve/e_l) - 1;
    fprintf("making Kp =  %0.3f\n", Kp_ll);
    Kp_target = Kp_ll / Kp_l;
    fprintf("ratio for lead lag:  %0.3f\n", Kp_target);
    z_ll = Kp_target/1000;
    p_ll = 1/1000;
    fprintf("to keep ratio= %0.3f, use zll = %0.03f  ,  pll= %0.03f\n", (z_ll/p_ll),z_ll, p_ll);
    %
    Gs_ll = (s + z_ll) / (s + p_ll);
    temp = Gs_ll * Gss_lead;
    tempp = temp/(1+temp);
    tf_temp = f_Fs_2_tf(tempp);
    h_poles = pole(tf_temp)
    h_zeros = zero(tf_temp)
    %
    %
        temp = (1417*(s+5)*(s+0.04759))/((s+2)*(s+4)*(s+6)*(s+8)*(s+0.001)*(s+15.18));
        tempp = temp/(1+temp);
        tf_temp = f_Fs_2_tf(tempp);
        h_poles = pole(tf_temp)
        h_zeros = zero(tf_temp)
    %

    temp_u = (355)/((s+2)*(s+4)*(s+6)*(s+8));
    u_tf = f_Fs_2_tf(temp_u);
    temp_c = (1417*(s+5)*(s+0.04759))/((s+2)*(s+4)*(s+6)*(s+8)*(s+0.001)*(s+14.679));
    c_tf = f_Fs_2_tf(temp_c);
    t_final = 50;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(u_tf, t_final, "b-");
    step(c_tf, t_final, "r-");
    grid on;
    legend({'Uncompensated', 'Compensated'},...
        Location='best', FontSize=14);
    templ = findall(gcf,'type','line');
    set([templ(5), templ(7)], 'linewidth',2);
    set([templ(1), templ(2)], 'linewidth',1);
    hold off;
%}

%{
    zeta = 0.5;
    theta = acos(zeta);
    os_given = f_zeta_2_os(zeta);
    fprintf("\ngiven zeta, theta=  %0.1f degree,  the %%OS=  %0.3f\n",...
        rad2deg(theta), 100*os_given);
   
    syms s;
    Gs = 1 / ((s + 2) * (s + 4) * (s + 6) * (s + 8));
    tf_initial = f_Fs_2_tf(Gs);
    f_rlocus(tf_initial); % use s = -1.53 + j2.65, K = 355
    K_u = 355;
    s_lag = -1.53 + 1j*2.65;
    KGs_u = K_u * Gs;
    temp = KGs_u / (1 + KGs_u);
    u_tf = f_Fs_2_tf(temp);

    Ts = 4 / abs(real(s_lag));
    fprintf("\nthe current Ts = %0.3f seconds\n", Ts);
    Ts_new = Ts - 0.5;
    fprintf("the design needs Ts = %0.3f\n", Ts_new);
    zeta_wn = -1 * 4 / Ts_new;
    fprintf("this means zeta * wn =  %0.3f\n", zeta_wn);
    s_op_i = -1 * tan(theta) * zeta_wn;
    s_op = zeta_wn + 1j * s_op_i;
    fprintf("operate now at s= %0.3f + j%0.3f\n",...
        real(s_op), imag(s_op));

    z_lead = 5;
    v_z = angle(s_op - (-1)*z_lead);
    v_1 = angle(s_op - (-1)*2);
    v_2 = angle(s_op - (-1)*4);
    v_3 = angle(s_op - (-1)*6);
    v_4 = angle(s_op - (-1)*8);
    anglez = v_z - (v_1 + v_2 + v_3 + v_4);
    fprintf("\nwithout the lead compensator pole, the angles are:  %0.3f degree\n", rad2deg(anglez));
    pc_lead_ang = pi + anglez;
    fprintf("the vector angle of the pole:  %0.3f degree\n", rad2deg(pc_lead_ang));
    p_lead = -1*real(s_op) + (imag(s_op)/tan(pc_lead_ang));
    fprintf("the poles is at:  %0.3f\n", -1*p_lead);

    G_lead = (s + z_lead) / (s + p_lead);
    Gs_l = Gs * G_lead;
    tf_lead = f_Fs_2_tf(Gs_l);
    f_rlocus(tf_lead);
    
    K_l = 1360;
    Gss_lead = K_l * Gss;
    Kp_l = limit(Gss_lead, s, 0);
    fprintf("Kp = %0.3f\n", Kp_l);
    e_l = 1 / (1 + Kp_l);
    fprintf("with lead error:  %0.3f\n", e_l);
    e_target = e_u/30;
    fprintf("make the error:  %0.3f\n", e_target);
    improve = e_l / e_target;
    fprintf("improve this much:  %0.3f\n", improve);
    e_ll = e_l / improve;
    fprintf("final error is:  %0.3f\n", e_ll);
    Kp_ll = (improve/e_l) - 1;
    fprintf("making Kp =  %0.3f\n", Kp_ll);
    Kp_target = Kp_ll / Kp_l;
    fprintf("ratio for lead lag:  %0.3f\n", Kp_target);
    z_ll = Kp_target/1000;
    p_ll = 1/1000;
    fprintf("to keep ratio= %0.3f, use zll = %0.03f  ,  pll= %0.03f\n", (z_ll/p_ll),z_ll, p_ll);

    t_final = 50;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(u_tf, t_final, "b-");
    grid on;
    %legend({'Uncompensated', 'Compensated'}, Location='best', FontSize=14);
    templ = findall(gcf,'type','line');
    %set([templ(5), templ(7)], 'linewidth',2);
    %set([templ(1), templ(2)], 'linewidth',1);
    hold off;



%}
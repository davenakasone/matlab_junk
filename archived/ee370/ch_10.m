%{
    chapter 10 : stability, nyquist

    1  :  class ex
    2  :  class ex
    3  :  class ex
    4  :  class ex
    5  :  class ex
    6  :  quiz, p1, i

    rlocus()
    rlocusplot()
    sgrid()
    step()
    Control System Designer
    nyquist()
    nicholas()

%}
format compact;
clear;
close all;
clc;
select = 6;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    z_1 = 1;
    z_2 = 1;
    k = 0.5;
    Ls = k / (s * (z_1*s + 1) * (z_2*s + 1));
    Ts = Ls / (Ls + 1);
    
    [tnum, tden] = numden(Ts); 
    T_num = sym2poly(tnum);
    T_den = sym2poly(tden);
    Tsys = tf(T_num, T_den)

    [lnum, lden] = numden(Ls); 
    L_num = sym2poly(lnum);
    L_den = sym2poly(lden);
    Lsys = tf(L_num, L_den)

    bode(Tsys);
    set(findall(gcf,'type','line'),'linewidth',2);
    f_rlocus(Lsys);

    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    %axis equal;
    nyquist(Tsys, 'r-');
    set(findall(gcf,'type','line'),'linewidth',2);
    set(groot,'defaultLineMarkerSize',15);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms alpha;
    syms s;
    syms z;

    G_c = (1 + z*s) / (1 + z*alpha*s);
    G_cc = subs(G_c, [alpha, z], [2, 0.2]);
    
    T_cc = G_cc;
    [t_num, t_den] = numden(T_cc);
    T_num = sym2poly(t_num);
    T_den = sym2poly(t_den);
    T_sys = tf(T_num, T_den);

    bode(T_sys);
    set(findall(gcf,'type','line'),'linewidth',2);
end

%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    zeta = 0.45; % min
    theta = acos(zeta);
    Kv = 20;
    K = 40;
    phi_m = deg2rad(45); % because phi_m = zeta/0.01

    Gs = 40 / (s * (s + 2));
    tf_init = f_Fs_2_tf(Gs);
    %f_rlocus(tf_init);
    Ts = Gs/(1+Gs);
    tf_f = f_Fs_2_tf(Ts);
    
    Gs_c = (s + 0.15) / (10*(s+0.15));
    tfg = f_Fs_2_tf(Gs_c);
    TT = Gs*Gs_c / (1 + Gs*Gs_c);
    tf_c = f_Fs_2_tf(TT);

    figure(1);
    bode(tf_init);
    set(findall(gcf,'type','line'),'linewidth',2);

    figure(2);
    axis equal;
    sgrid;
    nyquist(tf_f, 'r-', tf_c, 'b-');
    set(findall(gcf,'type','line'),'linewidth',2);

    figure(3);
    bode(tfg);
    set(findall(gcf,'type','line'),'linewidth',2);

    

    t_final = 10;
    figure(Position=[20, 20, 800, 800]);
    hold on;
    step(tf_f, t_final, "r--");
    step(tf_c, t_final, "b-");
    grid on;
    templ = findall(gcf,'type','line');
    set([templ(5), templ(7)], 'linewidth',3);
    hold off;

    pol_u = pole(tf_f)
    zer_u = zero(tf_f)
    pol_c = pole(tf_c)
    zer_c = zero(tf_c)

end

%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    K = 10;
    Gs = K / (s + 1)^2;
    tf_u = f_Fs_2_tf(Gs);

    %bode(tf_u);
    %nyquist(tf_u);

    syms w;
    w_c = 3;
    Gjw = subs(Gs, s, 1j*w_c);
    test = abs(Gjw); % should be abs==1
    test = 0-atan(w_c)-atan(w_c); %solve the angle
    p_m = pi + test; % can rotate this much before unstable
    fprintf("angle is:  %0.3f degree, can rotate:  %0.3f\n",...
        rad2deg(test), rad2deg(p_m));

    alpha = 2.5;
    tau = 6;
    Gs_lag = (1 + tau*s) / (1 + tau*alpha*s);
    tf_lag = f_Fs_2_tf(Gs_lag);
    Gs_c = Gs * Gs_lag;
    tf_c = f_Fs_2_tf(Gs_c);

    figure;
    hold on;
    bode(tf_u, "b--");
    bode(tf_lag, "r--");
    bode(tf_c, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("original", "lag", "compensated", location='best', fontsize=14);
    hold off;

    figure;
    hold on;
    axis equal;
    n = 1000; %// Define number of points on circle
    theta = linspace(0, 2*pi, n);
    x = cos(theta);
    y = sin(theta);
    nyquist(tf_u, "b--");
    nyquist(tf_lag, "r--");
    nyquist(tf_c, "g-");
    plot(x,y, "k:");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("original", "lag", "compensated", location='best', fontsize=14);
    hold off;

   test = atan(deg2rad(60));
   fprintf("rads:  %0.3f\n", test);

   test = 10/(1+ 1.7^2)
    
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    K = 10;
    Gs = K / (s + 1)^2;
    tf_u = f_Fs_2_tf(Gs);

    %bode(tf_u);
    %nyquist(tf_u);

    syms w;
    w_c = 3;
    Gjw = subs(Gs, s, 1j*w_c);
    test = abs(Gjw); % should be abs==1
    test = 0-atan(w_c)-atan(w_c); %solve the angle
    p_m = pi + test; % can rotate this much before unstable
    fprintf("angle is:  %0.3f degree, can rotate:  %0.3f\n",...
        rad2deg(test), rad2deg(p_m));

    alpha = 3;
    tau = 6;
    Gs_lead = (1 + alpha*tau*s) / (1 + tau*s);
    tf_lead = f_Fs_2_tf(Gs_lead);
    Gs_c = Gs * Gs_lead;
    tf_c = f_Fs_2_tf(Gs_c);

    figure;
    hold on;
    bode(tf_u, "b--");
    bode(tf_lead, "r--");
    bode(tf_c, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("original", "lead", "compensated", location='best', fontsize=14);
    hold off;

    figure;
    hold on;
    axis equal;
    n = 1000; %// Define number of points on circle
    theta = linspace(0, 2*pi, n);
    x = cos(theta);
    y = sin(theta);
    nyquist(tf_u, "b--");
    nyquist(tf_lead, "r--");
    nyquist(tf_c, "g-");
    plot(x,y, "k:");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("original", "lead", "compensated", location='best', fontsize=14);
    hold off;

   test = atan(deg2rad(60));
   fprintf("rads:  %0.3f\n", test);

   test = 10/(1+ 1.7^2)
    
end


%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    K_max = 250;
    K_min = 0.01;
    K = 100;
    Ls = 1 / (s * (s + 5)^2);
    KLs_max = K_max * Ls;
    KLs_min = K_min * Ls;
    KLs = K * Ls;
    tf_KL = f_Fs_2_tf(KLs);
    tf_KL_max = f_Fs_2_tf(KLs_max);
    tf_KL_min = f_Fs_2_tf(KLs_min);

    figure;
    hold on;
    bode(tf_KL_max, "b--");
    bode(tf_KL_min, "r--");
    bode(tf_KL, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("max K", "min K", "K = 100", location='best', fontsize=14);
    hold off;

    figure;
    hold on;
    axis equal;
    n = 100;
    theta = linspace(0, 2*pi, n);
    x = cos(theta);
    y = sin(theta);
    nyquist(tf_KL_max, "b--");
    nyquist(tf_KL_min, "r--");
    nyquist(tf_KL, "g-");
    plot(x,y, "k:");
    set(findall(gcf,'type','line'),'linewidth',2);
    legend("max K", "min K", "K = 100", location='best', fontsize=14);
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
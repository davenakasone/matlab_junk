%{
    chapter 10 : stability, nyquist

    1  :  class ex
    2  :  class ex
    3  :  class ex

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
select = 3;


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
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
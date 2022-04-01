%{
    chapter 8 : root locus

    1  :  ex8.1, basic complex vector
    2  :  sa8.1, use vector or sub
    3  :  sa8.2, poles, zeros, sub
    4  :  ex8.3, break away, break in  by differentiation + summings
    5  :  ex8.6, angles, + rlocus()
    6  :  hw_p9, sys1, design some systems, find values
    7  :  hw_p9, sys2, design some systems, find values
    8  :  hw_p12, get info
    9  :  hw_p14, basic root locus
    10 :  hw_p20, root locus
    11 :  hw_p21, root...locus
    12 :  hw_p25,
    13 :  hw_p28
    14 :  hw_p43
    15 :  class ex
    16 : class ex
    17 : continue # 16

    rlocus()
    rlocusplot()
    sgrid()
    step()
    Control System Designer

%}
format compact;
clear;
close all;
clc;
select = 17;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    Fs = (s + 1) / (s * (s + 2));
    val = -3 + 1j*4;
    
    vec1 = f_make_vec(-1, val);
    f_rec2pol(vec1);
    vec2 = f_make_vec(0, val);
    f_rec2pol(vec2);
    vec3 = f_make_vec(-2, val);
    f_rec2pol(vec3);
    
    fprintf("\nvalue of the point:\n");
    Fss = subs(Fs, s, val);
    f_rec2pol(Fss);
    fprintf("total of vectors:\n");
    M = abs(vec1) / (abs(vec2) * abs(vec3));
    tht = angle(vec1) - (angle(vec2) + angle(vec3));
    total = M * exp(1j*tht);
    f_rec2pol(total);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    ss = -7 + 1j * 9;
    Fs = ((s + 2) * (s + 4)) / (s * (s + 3) * (s + 6));
    subbing = subs(Fs, s, ss);
    fprintf("\n the substitution value:\n");
    f_rec2pol(subbing);
    
    z_1 = f_make_vec(-2, ss);
    z_2 = f_make_vec(-4, ss);
    p_1 = f_make_vec(0, ss);
    p_2 = f_make_vec(-3, ss);
    p_3 = f_make_vec(-6, ss);
    M = (abs(z_1) * abs(z_2)) / (abs(p_1) * abs(p_2) * abs(p_3));
    tht = (angle(z_1) + angle(z_2)) - (angle(p_1) + angle(p_2) + angle(p_3));
    fprintf("\n the vector value:\n");
    f_rec2pol(M * exp(1j*tht));
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
 
    ss = -3 + 1j * 0;
    Gs = (s + 2) / (s^2 + 4*s + 13);
    val = subs(Gs, s, ss);
    fprintf("\nvalue at point by subbing:\n");
    f_rec2pol(val);
    z_1 = f_make_vec(-2, ss);
    p_1 = f_make_vec(-2 + 1j * 3, ss);
    p_2 = f_make_vec(-2 - 1j * 3, ss);
    M = abs(z_1) / (abs(p_1) * abs(p_2));
    tht = angle(z_1) - angle(p_1) - angle(p_2);
    f_rec2pol(M * exp(1j*tht));
    
    K = 1 / M;
    KGs = subs(Gs, s, ss) * K;
    f_rec2pol(KGs);
    
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    GHs = ((s-3)*(s-5)) / ((s+1)*(s+2));
    K = 1 / GHs;
    K_d = diff(K, s);
    sigs = double(solve(K_d == 0, s))
    
    syms sig;
    eqn_zeros = (1/(sig-3)) + (1/(sig-5));
    eqn_poles = (1/(sig+1)) + (1/(sig+2));
    breaks = double(solve(eqn_zeros == eqn_poles, sig))
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    Gs = (s + 2) / ((s + 3) * (s^2 + 2*s + 2));
    %pretty(expand(Gs));
    my_tf = tf([1, 2], [1, 5, 8, 6]); % need to make an obj
    %rlocus(my_tf);
    
    rootz = roots([1, 5, 8, 6])
end


%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    syms K;
    Gs = ((s + 2) * (s + 1)) / ((s - 2) * (s - 1));
    totf = Gs / (1 + Gs);
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, 3, 2], [1, -3, 2]);
    character = (s - 2) * (s - 1) + K * (s + 2) * (s + 1);
    
    syms sig;
    eqn_zeros = (1/(sig + 2)) + (1/(sig + 1));
    eqn_poles = (1/(sig - 2)) + (1/(sig - 1));
    breaks = double(solve(eqn_zeros == eqn_poles, sig));
    sig_away = breaks(1);
    sig_in = breaks(2);
    K_away = solve(subs(character, s, sig_away)==0, K);
    K_in = solve(subs(character, s, sig_in)==0, K);
    fprintf("\na) break away:  %0.3f  ,  K = %0.3f\n", sig_away, K_away);
    fprintf("a) break in  :  %0.3f  ,  K = %0.3f\n\n", sig_in, K_in);
    
    roth = expand(subs(character, K, 1))
    KGs = K * Gs;
    Ts = KGs/(1 + KGs);
    %pretty(simplify(expand(Ts)));
    
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',20);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 7)
    syms s;
    syms K;
    Gs = ((s + 2) * (s + 1)) / (s^2 - 2*s + 2);
    totf = Gs / (1 + Gs);
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, 3, 2], [1, -2, 2]);
    character = (s^2 - 2*s + 2) + K * (s + 2) * (s + 1);
    
    syms sig;
    eqn_zeros = (1/(sig + 2)) + (1/(sig + 1));
    eqn_poles = (1/(sig + 1 + 1j*1)) + (1/(sig + 1 - 1j*1));
    breaks = double(solve(eqn_zeros == eqn_poles, sig))
    sig_away = breaks(1);
    sig_in = breaks(2);
    K_away = solve(subs(character, s, sig_away)==0, K);
    K_in = solve(subs(character, s, sig_in)==0, K);
    fprintf("\na) break away:  %0.3f  ,  K = %0.3f\n", sig_away, K_away);
    fprintf("a) break in  :  %0.3f  ,  K = %0.3f\n\n", sig_in, K_in);
    
    roth = expand(subs(character, K, 2/3))
    KGs = K * Gs;
    Ts = KGs/(1 + KGs);
    pretty(simplify(expand(Ts)));
    
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',20);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 8)
    syms s;
    syms sig;
    syms K;
    Gs = 1 / (s * (s + 5) * (s + 8));
    totf = Gs / (1 + Gs);
    pretty((expand(Gs)));
    my_tf = tf(1, [1, 13, 40, 0]);
    character = s^3 + 13*s^2 + 40*s + K;
    
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',20);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 9)
    syms s;
    syms sig
    syms K;
    Gs = (s + 3) / (s * (s + 1) * (s + 4) * (s + 6));
    KGs = K * Gs;
    Ts = Gs / (1 + Gs);
    KTs = KGs / (1 + KGs);
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, 3], [1, 11, 34, 24, 0]);
    pretty(simplify(expand(KTs)));
    character = s^4 + 11*s^3 + 34*s^2 + (K + 24)*s + 3*K;
    
    brk_eqn = -1/Gs;
    brk_eqnn = diff(brk_eqn, s);
    brk_pt = double(solve(brk_eqnn == 0, s))
    KK = double(subs(brk_eqn, s, brk_pt(1)))
    
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 10)
    syms s;
    syms sig;
    syms K;
    Gs = (s^2 - 2*s + 2) / ((s + 1) * (s + 3) * (s + 4) * (s + 5));
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, -2, 2], [1, 13, 59, 107, 60]);
    KGs = K * Gs;
    Ts = Gs / (1 + Gs);
    KTs = KGs / (1 + KGs);
    pretty(simplify(expand(Ts)));
    pretty(simplify(expand(KTs)));
    
    brk_eqn = -1/Gs;
    brk_eqnn = diff(brk_eqn, s);
    brk_pt = double(solve(brk_eqnn == 0, s))
    KKa = double(subs(brk_eqn, s, brk_pt(1)))
    KKb = double(subs(brk_eqn, s, brk_pt(5)))
    
    zeta = -1 * log(0.2) / sqrt(pi^2 + (log(0.2))^2)
    
    character = s^4 + 13*s^3 + (K + 59)*s^2 + (107 - 2*K)*s + 2*K + 60;
    eqn = subs(character, K, 13);
    rootz = roots([1, 13, 72, 81, 86]);
    
    Tt = subs(KTs, K, 13);
    pretty(simplify(Tt));
    my_g = tf([13, -26, 26], [1, 13, 72, 81, 86]);
    step(my_g);
    set(findall(gcf,'type','line'),'linewidth',2);
    
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 11)
    syms K;
    syms s;
    Gs = ((s + 2) * (s + 3)) / (s * (s + 1));
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, 5, 6], [1, 1, 0]);
    KGs = K * Gs;
    Ts = Gs / (1 + Gs);
    pretty(simplify(expand(Ts)));
    KTs = KGs / (1 + KGs);
    pretty(simplify(expand(KTs)));
    character = (1 + K)*s^2 + (5*K + 1)*s + 6*K;
    
    zeta_fun = (1/2)*(1/sqrt((6*K)/(K+1)))*((5*K+1)/(K+1));
    d_zeta_fun = diff(zeta_fun, K);
    minz = solve(d_zeta_fun==0, K);
    w_n = double(sqrt(6*minz/(minz+1)))
    zeta = double(subs(zeta_fun, K, minz));
    o_s = exp((-zeta*pi)/sqrt(1-zeta^2))*100
    t_settle = 4 / (zeta * w_n)
    t_peak = pi/(w_n * sqrt(1-zeta^2))
    
    fun_T = subs(KTs, K, 1/3);
    pretty(simplify(expand(fun_T)));
    fun_G = subs(KGs, K, 1/3);
    pretty(simplify(expand(fun_G)));
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 12)
    syms s;
    syms K;
    Gs = (s^2 + 4*s +5) / ((s^2 + 2*s +5) * (s + 3) * (s + 4));
    pretty(simplify(expand(Gs)));
    my_tf = tf([1, 4, 5], [1, 9, 31, 59, 60]);
    
    w_nn = sqrt((5*K+5)/(K+1));
    zetaa = (1/2)*((4*K+2)/(K+1))*(1/w_nn);
    k_k = double(solve(w_nn*sqrt(1-zetaa^2)==pi, K))
    wn_a = subs(w_nn, K, k_k(1));
    zeta_a = subs(zetaa, K, k_k(1));
    wn_b = subs(w_nn, K, k_k(2));
    zeta_b = subs(zetaa, K, k_k(2));
    tp_a = double(pi/(wn_a*sqrt(1-zeta_a^2)))
    tp_b = double(pi/(wn_b*sqrt(1-zeta_b^2)))
    
    KGs = 10.8 * Gs;
    KTs = KGs / (1 + KGs);
    pretty(simplify(expand(KTs)));
    my_t = tf(54.*[1, 4, 5], [5, 45, 209, 511, 570]);
    step(my_t);
    set(findall(gcf,'type','line'),'linewidth',2);
    
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 13)
    syms K;
    syms sig;
    syms s;
    Gs = 1 / ((s + 2) * (s + 3));
    Hs = (s^2 - 4*s + 8) / (s^2 + 2*s + 5);
    GHs = Gs * Hs;
    pretty(simplify(expand(GHs)));
    my_tf = tf([1, -4, 8], [1, 7, 21, 37, 30]);
    
    my_zeros = roots([1, -4, 8]);
    z_1 = my_zeros(1);
    z_2 = my_zeros(2);
    my_poles = roots([1, 7, 21, 37, 30]);
    p_1 = my_poles(1);
    p_2 = my_poles(2);
    p_3 = my_poles(3);
    p_4 = my_poles(4);
    
    z2_a1 = angle(f_make_vec(z_2, z_1));
    p1_a1 = angle(f_make_vec(p_1, z_1));
    p2_a1 = angle(f_make_vec(p_2, z_1));
    p3_a1 = angle(f_make_vec(p_3, z_1));
    p4_a1 = angle(f_make_vec(p_4, z_1));
    z1_a1 = rad2deg(pi + p1_a1 + p2_a1 + p3_a1 + p4_a1 - z2_a1)
    
    z1_a2 = angle(f_make_vec(z_1, z_2));
    p1_a2 = angle(f_make_vec(p_1, z_2));
    p2_a2 = angle(f_make_vec(p_2, z_2));
    p3_a2 = angle(f_make_vec(p_3, z_2));
    p4_a2 = angle(f_make_vec(p_4, z_2));
    z2_a2 = rad2deg(pi + p1_a2 + p2_a2 + p3_a2 + p4_a2 - z1_a2)
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 14)
    syms K;
    syms s;
    Gs = (61.73/(s^2 + 11.11*s + 61.73)) * (1/(s+10)^3);
    pretty(expand(simplify(Gs)));
    my_tf = tf(6173, [100, 4111, 69503, 618490, 2962900, 6173000]);
    KTs = (K*Gs)/(1 + K*Gs);
    pretty(expand(simplify(KTs)));
    
    zeta = -1 * log(0.2) / sqrt(pi^2 + (log(0.2))^2)
    
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 15)
    syms s;
    syms K;
    zeta = -1 * log(0.16) / sqrt(pi^2 + (log(0.16))^2)%0.5039
    
    Gs = 1/(s * (s + 4) * (s + 6));
    pretty(expand(Gs));
    my_tf = tf(1, [1, 10, 24, 0]);
    
    Ts = Gs / (1 + Gs);
    KGs = K * Gs;
    KTs = KGs / (1 + KGs);
    
    myT = subs(KTs, K, 43.35)
    
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 16)
    syms s;
    syms K;
    zeta = -1 * log(0.3) / sqrt(pi^2 + (log(0.3))^2)%
    fprintf("\nzeta = %0.3f\n", zeta);
    
    Gs = 1/ (s * (s + 4) * (s + 5));
    %pretty(expand(Gs));
    my_tf = tf(1, [1, 9, 20, 0]);
    breaks = double(solve(diff(1/Gs)==0, s)); 
    fprintf("\nbreak away is:  %0.3f\n", breaks(2));
    T_s = 4 / (zeta * abs(-0.93 + 1j * 2.24));
    fprintf("\nw_n = %0.3f\n", T_s);
    fprintf("\nuncompensated settle time:  %0.3f\n", T_s);
    T_s = 4 / (0.93);
    fprintf("\nuncompensated settle time:  %0.3f\n", T_s);
    
    %{
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(my_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 17)
    zeta = -1 * log(0.3) / sqrt(pi^2 + (log(0.3))^2)
    zeta_wn = 2 * 0.93;  % must be twice to reduce settling time by 2
    known = -2.014 + 1j*5.252;
    syms pc;
    p1 = f_make_vec(pc, known);
    p2 = f_make_vec(-6, known);
    p3 = f_make_vec(-6, known);
    p4 = f_make_vec(-6, known);
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>

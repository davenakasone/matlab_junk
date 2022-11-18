%{
    final review

    1  :  Routh, special case 1, 0 in a column, use epsilon...try preparation also
    2  :  Routh, special case 2, row becomes 0
    3  :  sketching root locus
    4  :  root locus design, probably a lead
    5  :  lag, nyquist

%}
format compact;
clear;
close all;
clc;
select = 3;


%------------------------------------------------------------------------------------------
if (select == 1)
    % given char: s^5 + 2*s^4 + 3*s^3 + 6*s^2 + 5*s + 3    CHANGE
    syms s;
    syms ep;
    characteristic = (s^5 + 2*s^4 + 3*s^3 + 6*s^2 +5*s + 3);
    characteristic_p = sym2poly(characteristic)
    quick_check = roots(characteristic_p)

    a_5 = 1; a_3 = 3; a_1 = 5;   %  0
    a_4 = 2; a_2 = 6; a_0 = 3;   %  0
    % c_11     c_12    % 0       %  0
    % c_21     c_22    % 0
    % c_31     c_32    % 0

    % first row, s^3
    c_11 = f_routh(a_5, a_3, a_4, a_2, a_4) % it == 0, so make it epsilon
    c_11 = ep;
    c_12 = f_routh(a_5, a_1, a_4, a_0, a_4);

    % second row, s^2
    c_21 = f_routh(a_4, a_2, c_11, c_12, c_11);
    c_22 = f_routh(a_4, a_0, c_11, 0, c_11);

    % third row, s^1
    c_31 = f_routh(c_11, c_12, c_21, c_22, c_21);
    c_32 = f_routh(c_11, 0, c_21, 0, c_21);

    % fourth row
    c_41 = f_routh(c_21, c_22, c_31, c_32, c_31);

    epz = linspace(-5, 5, 100);
    figure;
    grid on;
    hold on;
    plot(epz, subs(c_21, ep, epz), "r-", LineWidth=2);
    plot(epz, subs(c_31, ep, epz), "b-", LineWidth=2);
    hold off;

end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    Cs = s^5 + 7*s^4 + 6*s^3 + 42*s^2 + 8*s + 56;
    check_roots = roots(sym2poly(Cs))
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    K = 1;
    Gs = K / (s * (s + 20) * (s+40));
    tf_Gs = f_Fs_2_tf(Gs);
    f_rlocus(tf_Gs);

    syms sig;
    KK = -1/subs(Gs, s, sig);
    dKK = diff(KK, sig);
    siga = double(solve(dKK==0, sig))
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    Gs = 1 / (s * (s + 1));
    zeta = 1/sqrt(2);
    w_n = 3*sqrt(2);
    domA = -1 * zeta * w_n + w_n * sqrt(zeta^2 - 1);
    domB = -1 * zeta * w_n - w_n * sqrt(zeta^2 - 1);
    fprintf("so design at s=  %0.3f +/- j%0.3f\n", real(domA), imag(domA));
    %{
    v_1 = angle(domA - 0);
    v_2 = angle(domA +10);
    v_3 = angle(domA +50);
    v_4 = angle(domA +120);
    theta = pi - [v_1 + v_2 + v_3 + v_4];
    fprintf("sum of these poles means zero angle is:  %0.3f deg\n", rad2deg(theta));
    z_c = (imag(domA) + tan(theta)*real(domA)) / tan(theta);
    fprintf("this means zc=  %0.3f\n", z_c);
    %}
    Gs_lead = (s +3)/(s + 18.0135);
    Ls = Gs * Gs_lead;
    tf_Ls = f_Fs_2_tf(Ls);
    f_rlocus(tf_Ls);

    comp_K = 1/abs(subs(Ls, s, domA));
    fprintf("K=%0.3f\n", comp_K);
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    Ls = 475/ (s+5)^2;
    Ls = Ls * ((1+s*0.58)/(1+1.97*s))
    tf_Ls = f_Fs_2_tf(Ls);
    bode(tf_Ls)

    w = tan(deg2rad(65))*5
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
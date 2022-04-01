%{
    1: exp1
    2: exp2
    3: exp3
%}
format compact;
close all;
clear;
clc;
select = 3;


%-------------------------------------------------------------------------------------
if select == 1
    V_s = 124.9;
    S_abs = 0.003;
    P_s = 0.003;
    
    I_s = S_abs / V_s
    Y = I_s / V_s
    Q_s = -1 * Y * V_s^2
end


%-------------------------------------------------------------------------------------
if select == 2
    Y = 1.9231e-07;
    V_s = 124.8;
    I_s = 2.055;
    S_s_abs = 256.6;
    P_s = 38.27;
    Q_s = 252.9;
    I_r = 2.055;
    
    R_line = P_s / I_r^2
    X_line = (Q_s + (Y/2)*V_s^2) / I_r^2
    Z_line = R_line + 1j * X_line
    [A, B, C, D] = f_line_med_ABCD(Z_line, Y)
    
    check_Vs = Z_line * I_r;
    f_mdri("Vs", check_Vs, 1);
    check_Is = D * I_r;
    f_mdri("Is", check_Is, 1);
    check_Vr = V_s / A;
    f_mdri("Vr", check_Vr, 1);
end


%-------------------------------------------------------------------------------------
if select == 3
    %[A, B, C, D] = f_line_med_ABCD(z_line, y_line);
    z_line = 9.0622 + 1j * 59.8863;
    y_line = 1.9231e-7;
    y2z2 = (2 / y_line) / (1j);
    arr_siz = 21;
    V_s = 124.6;
    r_load = [ 1200, 600, 400, 300, 240, 200, 171, 150, 133, 120, 109, 100, 92, 86, 80, 75, 71, 67, 63, 60, 57];
    I_s = zeros(1, arr_siz);
    P_s = zeros(1, arr_siz);
    Q_s = zeros(1, arr_siz);
    I_r = zeros(1, arr_siz);
    V_r = zeros(1, arr_siz); 
    P_r = zeros(1, arr_siz);
    eff_eta = zeros(1, arr_siz);
    volt_reg = zeros(1, arr_siz);
    temp_Pr_vr = 0;
    P_r_vr = 0;
    temp_Pr_eff = 0;
    P_r_eff = 0;
    fprintf("\nR_load  | V_s (V) | I_s (A) | P_s (W) | Q_s (VAR) | I_r (A)  | V_r (V) | P_r (W) | eff (%%) | VR (%%)\n");
    fprintf("-----------------------------------------------------------------------------------------------------\n");
    %
    for ii = 1:1:arr_siz
        A = [ y2z2 , -y2z2          , 0                   ;
              -y2z2, 2*y2z2 + z_line, -y2z2               ;
              0    , -y2z2          , y2z2 + r_load(1, ii)];
        Y = [ V_s; 0; 0];
        I = A \ Y;
        I_s(1, ii) = I(1,1);
        I_r(1, ii) = I(3,1);
        V_r(1, ii) = abs(I_r(1, ii)) * r_load(1, ii);
        P_s(1, ii) = real(V_s * conj(I_s(1, ii)));
        Q_s(1, ii) = imag(V_s * conj(I_s(1, ii)));
        P_r(1, ii) = real(I_r(1, ii))^2 * r_load(1, ii);
        eff_eta(1, ii) = (P_r(1, ii) / P_s(1, ii) ) * 100;
        volt_reg(1, ii) = ((V_s - V_r(1, ii)) / V_r(1, ii)) * 100;
        fprintf(" %6.1f | %7.2f | %7.2f | %7.2f | %9.2f | %7.2f  | %7.2f | %7.2f | %7.2f | %7.2f \n",...
            r_load(1, ii), V_s, I_s(1, ii), P_s(1, ii), Q_s(1, ii), I_r(1, ii), V_r(1, ii), P_r(1, ii), eff_eta(1, ii), volt_reg(1,ii));
        if temp_Pr_eff < P_r(1, ii)
            temp_Pr_eff = P_r(1, ii);
            P_r_eff = ii;
        end
        if temp_Pr_vr < P_r(1, ii)
            temp_Pr_vr = P_r(1, ii);
            P_r_vr = ii;
        end
    end
    %
    %
    figure('Position', [20, 20, 800, 800]);
    grid on;
    hold on;
    xlabel("power to load, P_r (W)", 'fontsize', 12);
    ylabel("percent (%)", 'fontsize', 12);
    title("voltage regulation and efficiency vs load power", 'fontsize', 18);
    plot(P_r(1,:), volt_reg(1,:), 'g-', 'linewidth', 2);
    plot(P_r(1,:), eff_eta(1,:), 'r-', 'linewidth', 2);
    plot(P_r(P_r_eff), eff_eta(1, P_r_eff), 'rx', 'markersize', 10, 'linewidth', 3);
    plot(P_r(P_r_vr), volt_reg(1, P_r_vr), 'gx', 'markersize', 10, 'linewidth', 3);
    legend('VR', 'eff','opt eff', 'opt vr' , 'fontsize', 12);
    a = P_r(1, P_r_vr)
    b = P_r_eff
    %
end
%%%%%%%%~~~~~~~~~END> 
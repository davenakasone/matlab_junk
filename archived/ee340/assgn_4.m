%{
    1  :  basic test, make the interval
    2  :  basic test, make 5% noise
    3  :  get the basic info
    4  :  variation applied
%}
close all;
clear all;
clc;
select = 4;


%-------------------------------------------------------------------------------------
if select == 1
    %{
    onez = 1:1:9;
    tenz = 10:10:990;
    thousandz = 1000:1000:20000;
    temp = cat(2, onez, tenz);
    inpz = cat(2, temp, thousandz)
    size(inpz)
    %}
    inpz = [1:1:9, 10:10:990, 1000:1000:20000]
    sizz = size(inpz)
    sizzz = sizz(1,2)
end

%-------------------------------------------------------------------------------------
if select == 2
    volts = 100;
    five_percent = 0.05 * volts;
    noise = zeros(1, 128);
    for ii = 1:1:128
        if mod(ii, 2)
            noise(1,ii) = five_percent * rand(1);
        else
            noise(1,ii) = -1 * five_percent * rand(1);
        end
    end
    temp = noise + volts
end


%-------------------------------------------------------------------------------------
if select == 3
    freq_v = 60; % Hz
    w = 2 * pi * freq_v; % rad / sec
    V_s = 345e3 / sqrt(3); % V
    S_abs = 300e6;
    theta_a = deg2rad(0);
    theta_b = deg2rad(30);
    theta_c = deg2rad(-30);
    Z_line = 200 * exp(1j * deg2rad(80)) % omhs
    Y_shunt = 1j * 0.001; % 1 / ohms
    [A, B, C, D] = f_line_med_ABCD(Z_line, Y_shunt);
    
    V_th = V_s / A;
    f_mdri("V_th", V_th, 1);
    %Z_th = f_para(Z_line, 2/Y_shunt);
    Z_th = 28.5 + 1j*179.5
    f_mdri("Z_th", Z_th, 1);
    
    z_load_absz = [1:1:9, 10:10:990, 1000:1000:20000];
    sizz = size(z_load_absz);
    size_z = sizz(1,2);
    
    z_load_a = zeros(1, size_z);
    z_load_b = zeros(1, size_z);
    z_load_c = zeros(1, size_z);
    V_r_a = zeros(1, size_z);
    V_r_b = zeros(1, size_z);
    V_r_c = zeros(1, size_z);
    P_r_a = zeros(1, size_z);
    P_r_b = zeros(1, size_z);
    P_r_c = zeros(1, size_z);
    V_r_a_max = 0;
    V_r_b_max = 0;
    V_r_c_max = 0;
    P_r_a_max = 0;
    P_r_b_max = 0;
    P_r_c_max = 0;
    z_load_a_max = 0;
    z_load_b_max = 0;
    z_load_c_max = 0;
    
    for ii = 1:1:size_z
        z_load_a(1, ii) = z_load_absz(ii) * exp(1j * theta_a);
        z_load_b(1, ii) = z_load_absz(ii) * exp(1j * theta_b);
        z_load_c(1, ii) = z_load_absz(ii) * exp(1j * theta_c);
        temp_z_a = abs(z_load_a(1, ii)) / abs(z_load_a(1, ii) + Z_th);
        temp_z_b = abs(z_load_b(1, ii)) / abs(z_load_b(1, ii) + Z_th);
        temp_z_c = abs(z_load_c(1, ii)) / abs(z_load_c(1, ii) + Z_th);
        V_r_a(1, ii) = (abs(V_th) * temp_z_a);
        V_r_b(1, ii) = (abs(V_th) * temp_z_b);
        V_r_c(1, ii) = (abs(V_th) * temp_z_c);
        
        temp_p_a = V_r_a(1, ii)^2 / abs(z_load_a(1, ii));
        temp_p_b = V_r_b(1, ii)^2 / abs(z_load_b(1, ii));
        temp_p_c = V_r_c(1, ii)^2 / abs(z_load_c(1, ii));
        P_r_a(1, ii) = temp_p_a * cos(theta_a);
        P_r_b(1, ii) = temp_p_b * cos(theta_b);
        P_r_c(1, ii) = temp_p_c * cos(theta_c);
        
        if P_r_a(1, ii) > P_r_a_max
            V_r_a_max = V_r_a(1, ii);
            P_r_a_max = P_r_a(1, ii);
            z_load_a_max = z_load_absz(ii);
        end
        if P_r_b(1, ii) > P_r_b_max
            V_r_b_max = V_r_b(1, ii);
            P_r_b_max = P_r_b(1, ii);
            z_load_b_max = z_load_absz(ii);
        end
        if P_r_c(1, ii) > P_r_c_max
            V_r_c_max = V_r_c(1, ii);
            P_r_c_max = P_r_c(1, ii);
            z_load_c_max = z_load_absz(ii);
        end
    end
    %
    figure('Position', [20, 20, 800, 800])
    hold on;
    grid on;
    title('Power vs Voltage, received', 'FontSize', 20);
    xlabel('Power received in MW', 'fontsize', 16);
    ylabel('Voltage received in kV', 'fontsize', 16);
    plot(P_r_a ./ 1e6, V_r_a ./ 1e3, 'r-', 'linewidth', 2);
    plot(P_r_b ./ 1e6, V_r_b ./ 1e3, 'b-', 'linewidth', 2);
    plot(P_r_c ./ 1e6, V_r_c ./ 1e3, 'g-', 'linewidth', 2);
    plot(P_r_a_max / 1e6, V_r_a_max / 1e3, 'ro', 'markersize', 10, 'linewidth', 3);
    plot(P_r_b_max / 1e6, V_r_b_max / 1e3, 'bo', 'markersize', 10, 'linewidth', 3);
    plot(P_r_c_max / 1e6, V_r_c_max / 1e3, 'go', 'markersize', 10, 'linewidth', 3);
    legend('angle(0)', 'angle(30)', 'angle(-30)',...
        'max', 'max', 'max', 'fontsize', 16, 'Location', 'best');
    text_a = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_a_max, P_r_a_max/1e6, V_r_a_max/1e3);
    t_a = text(100, 20, text_a); t_a.Color = 'red'; t_a.FontSize = 14;
    
    text_b = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_b_max, P_r_b_max/1e6, V_r_b_max/1e3);
    t_b = text(40, 20, text_b); t_b.Color = 'blue'; t_b.FontSize = 14;
    
    text_c = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_c_max, P_r_c_max/1e6, V_r_c_max/1e3);
    t_c = text(140, 50, text_c); t_c.Color = 'green'; t_c.FontSize = 14;
    %
    fprintf("\n\n\t\tPowers in (MW), Voltages in (kV)\n");
    fprintf("|Z_load|  :     V_r_0  :     P_r_0   :     V_r_30  :     P_r_30  :    V_r_-30  :   P_r_-30\n");
    fprintf("-------------------------------------------------------------------------------------------\n");
    for ii = 1:1:size_z
        fprintf("%6d    :  %9.3f :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f\n",...
            z_load_absz(ii),...
            V_r_a(1, ii) / 1e3, P_r_a(ii) / 1e6,...
            V_r_b(1, ii) / 1e3, P_r_b(ii) / 1e6,...
            V_r_c(1, ii) / 1e3, P_r_c(ii) / 1e6...
            );
    end
end


%-------------------------------------------------------------------------------------
if select == 4
    freq_v = 60; % Hz
    w = 2 * pi * freq_v; % rad / sec
    V_s = 345e3 / sqrt(3); % V
    S_abs = 300e6;
    theta_a = deg2rad(0);
    theta_b = deg2rad(30);
    theta_c = deg2rad(-30);
    Z_line = 200 * exp(1j * deg2rad(80)) % omhs
    Y_shunt = 1j * 0.001; % 1 / ohms
    [A, B, C, D] = f_line_med_ABCD(Z_line, Y_shunt);
    
    V_th = V_s / A;
    f_mdri("V_th", V_th, 1);
    %Z_th = f_para(Z_line, 2/Y_shunt);
    Z_th = 28.5 + 1j*179.5
    f_mdri("Z_th", Z_th, 1);
    
    z_load_absz = [1:1:9, 10:10:990, 1000:1000:20000];
    sizz = size(z_load_absz);
    size_z = sizz(1,2);
    
    z_load_a = zeros(1, size_z);
    z_load_b = zeros(1, size_z);
    z_load_c = zeros(1, size_z);
    V_r_a = zeros(1, size_z);
    V_r_a_down5 = zeros(1, size_z);
    V_r_a_up5 = zeros(1, size_z);
    V_r_b = zeros(1, size_z);
    V_r_b_down5 = zeros(1, size_z);
    V_r_b_up5 = zeros(1, size_z);
    V_r_c = zeros(1, size_z);
    V_r_c_down5 = zeros(1, size_z);
    V_r_c_up5 = zeros(1, size_z);
    P_r_a = zeros(1, size_z);
    P_r_a_down5 = zeros(1, size_z);
    P_r_a_up5 = zeros(1, size_z);
    P_r_b = zeros(1, size_z);
    P_r_b_down5 = zeros(1, size_z);
    P_r_b_up5 = zeros(1, size_z);
    P_r_c = zeros(1, size_z);
    P_r_c_down5 = zeros(1, size_z);
    P_r_c_up5 = zeros(1, size_z);
    V_r_a_max = 0;
    V_r_a_max_up5 = 0;
    V_r_a_max_down5 = 0;
    V_r_b_max = 0;
    V_r_b_max_up5 = 0;
    V_r_b_max_down5 = 0;
    V_r_c_max = 0;
    V_r_c_max_up5 = 0;
    V_r_c_max_down5 = 0;
    P_r_a_max = 0;
    P_r_a_max_up5 = 0;
    P_r_a_max_down5 = 0;
    P_r_b_max = 0;
    P_r_b_max_up5 = 0;
    P_r_b_max_down5 = 0;
    P_r_c_max = 0;
    P_r_c_max_up5 = 0;
    P_r_c_max_down5 = 0;
    z_load_a_max = 0;
    z_load_b_max = 0;
    z_load_c_max = 0;
    
    for ii = 1:1:size_z
        z_load_a(1, ii) = z_load_absz(ii) * exp(1j * theta_a);
        z_load_b(1, ii) = z_load_absz(ii) * exp(1j * theta_b);
        z_load_c(1, ii) = z_load_absz(ii) * exp(1j * theta_c);
        temp_z_a = abs(z_load_a(1, ii)) / abs(z_load_a(1, ii) + Z_th);
        temp_z_b = abs(z_load_b(1, ii)) / abs(z_load_b(1, ii) + Z_th);
        temp_z_c = abs(z_load_c(1, ii)) / abs(z_load_c(1, ii) + Z_th);
        V_r_a(1, ii) = (abs(V_th) * temp_z_a);
        V_r_a_down5(1, ii) = 0.95 * V_r_a(1, ii);
        V_r_a_up5(1, ii) = 1.05 * V_r_a(1, ii);
        V_r_b(1, ii) = (abs(V_th) * temp_z_b);
        V_r_b_down5(1, ii) = 0.95 * V_r_b(1, ii);
        V_r_b_up5(1, ii) = 1.05 * V_r_b(1, ii);
        V_r_c(1, ii) = (abs(V_th) * temp_z_c);
        V_r_c_down5(1, ii) = 0.95 * V_r_c(1, ii);
        V_r_c_up5(1, ii) = 1.05 * V_r_c(1, ii);
        
        temp_p_a = V_r_a(1, ii)^2 / abs(z_load_a(1, ii));
        temp_p_a_down5 = V_r_a_down5(1, ii)^2 / abs(z_load_a(1, ii));
        temp_p_a_up5 = V_r_a_up5(1, ii)^2 / abs(z_load_a(1, ii));
        
        temp_p_b = V_r_b(1, ii)^2 / abs(z_load_b(1, ii));
        temp_p_b_down5 = V_r_b_down5(1, ii)^2 / abs(z_load_b(1, ii));
        temp_p_b_up5 = V_r_b_up5(1, ii)^2 / abs(z_load_b(1, ii));
        
        temp_p_c = V_r_c(1, ii)^2 / abs(z_load_c(1, ii));
        temp_p_c_down5 = V_r_c_down5(1, ii)^2 / abs(z_load_c(1, ii));
        temp_p_c_up5 = V_r_c_up5(1, ii)^2 / abs(z_load_c(1, ii));
        
        P_r_a(1, ii) = temp_p_a * cos(theta_a);
        P_r_a_down5(1, ii) = temp_p_a_down5 * cos(theta_a);
        P_r_a_up5(1, ii) = temp_p_a_up5 * cos(theta_a);
        
        P_r_b(1, ii) = temp_p_b * cos(theta_b);
        P_r_b_down5(1, ii) = temp_p_b_down5 * cos(theta_b);
        P_r_b_up5(1, ii) = temp_p_b_up5 * cos(theta_b);
        
        P_r_c(1, ii) = temp_p_c * cos(theta_c);
        P_r_c_down5(1, ii) = temp_p_c_down5 * cos(theta_c);
        P_r_c_up5(1, ii) = temp_p_c_up5 * cos(theta_c);
        
        if P_r_a(1, ii) > P_r_a_max
            V_r_a_max = V_r_a(1, ii);
            P_r_a_max = P_r_a(1, ii);
            z_load_a_max = z_load_absz(ii);
        end
        if P_r_a_down5(1, ii) > P_r_a_max_down5
            V_r_a_max_down5 = V_r_a_down5(1, ii);
            P_r_a_max_down5 = P_r_a_down5(1, ii);
        end
        if P_r_a_up5(1, ii) > P_r_a_max_up5
            V_r_a_max_up5 = V_r_a_up5(1, ii);
            P_r_a_max_up5 = P_r_a_up5(1, ii);
        end
        
        if P_r_b(1, ii) > P_r_b_max
            V_r_b_max = V_r_b(1, ii);
            P_r_b_max = P_r_b(1, ii);
            z_load_b_max = z_load_absz(ii);
        end
        if P_r_b_down5(1, ii) > P_r_b_max_down5
            V_r_b_max_down5 = V_r_b_down5(1, ii);
            P_r_b_max_down5 = P_r_b_down5(1, ii);
        end
        if P_r_b_up5(1, ii) > P_r_b_max_up5
            V_r_b_max_up5 = V_r_b_up5(1, ii);
            P_r_b_max_up5 = P_r_b_up5(1, ii);
        end
        
        if P_r_c(1, ii) > P_r_c_max
            V_r_c_max = V_r_c(1, ii);
            P_r_c_max = P_r_c(1, ii);
            z_load_c_max = z_load_absz(ii);
        end
        if P_r_c_down5(1, ii) > P_r_c_max_down5
            V_r_c_max_down5 = V_r_c_down5(1, ii);
            P_r_c_max_down5 = P_r_c_down5(1, ii);
        end
        if P_r_c_up5(1, ii) > P_r_c_max_up5
            V_r_c_max_up5 = V_r_c_up5(1, ii);
            P_r_c_max_up5 = P_r_c_up5(1, ii);
        end
    end
    %
    figure('Position', [20, 20, 800, 800])
    hold on;
    grid on;
    axis padded;
    title('Power vs Voltage, received  with +/- 5%', 'FontSize', 20);
    xlabel('Power received in MW', 'fontsize', 16);
    ylabel('Voltage received in kV', 'fontsize', 16);
    plot(P_r_a ./ 1e6, V_r_a ./ 1e3, 'r-', 'linewidth', 2);
    plot(P_r_a_down5 ./ 1e6, V_r_a_down5 ./ 1e3, 'r:', 'linewidth', 2);
    plot(P_r_a_up5 ./ 1e6, V_r_a_up5 ./ 1e3, 'r--', 'linewidth', 2);
    plot(P_r_b ./ 1e6, V_r_b ./ 1e3, 'b-', 'linewidth', 2);
    plot(P_r_b_down5 ./ 1e6, V_r_b_down5 ./ 1e3, 'b:', 'linewidth', 2);
    plot(P_r_b_up5 ./ 1e6, V_r_b_up5 ./ 1e3, 'b--', 'linewidth', 2);
    plot(P_r_c ./ 1e6, V_r_c ./ 1e3, 'g-', 'linewidth', 2);
    plot(P_r_c_down5 ./ 1e6, V_r_c_down5 ./ 1e3, 'g:', 'linewidth', 2);
    plot(P_r_c_up5 ./ 1e6, V_r_c_up5 ./ 1e3, 'g--', 'linewidth', 2);
    
    plot(P_r_a_max / 1e6, V_r_a_max / 1e3, 'ro', 'markersize', 10, 'linewidth', 3);
    plot(P_r_a_max_down5 / 1e6, V_r_a_max_down5 / 1e3, 'ro', 'markersize', 10, 'linewidth', 3);
    plot(P_r_a_max_up5 / 1e6, V_r_a_max_up5 / 1e3, 'ro', 'markersize', 10, 'linewidth', 3);
    
    plot(P_r_b_max / 1e6, V_r_b_max / 1e3, 'bo', 'markersize', 10, 'linewidth', 3);
    plot(P_r_b_max_down5 / 1e6, V_r_b_max_down5 / 1e3, 'bo', 'markersize', 10, 'linewidth', 3);
    plot(P_r_b_max_up5 / 1e6, V_r_b_max_up5 / 1e3, 'bo', 'markersize', 10, 'linewidth', 3);
    
    plot(P_r_c_max / 1e6, V_r_c_max / 1e3, 'go', 'markersize', 10, 'linewidth', 3);
    plot(P_r_c_max_down5 / 1e6, V_r_c_max_down5 / 1e3, 'go', 'markersize', 10, 'linewidth', 3);
    plot(P_r_c_max_up5 / 1e6, V_r_c_max_up5 / 1e3, 'go', 'markersize', 10, 'linewidth', 3);
    
    legend('angle(0)', 'angle(0) - 5%', 'angle(0) + 5%',...
        'angle(30)', 'angle(30) - 5%', 'angle(30) + 5%',...
        'angle(-30)', 'angle(-30) - 5%', 'angle(-30) + 5%',...
        'fontsize', 16, 'Location', 'best');
    %
    text_a = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_a_max, P_r_a_max/1e6, V_r_a_max/1e3);
    t_a = text(100, 20, text_a); t_a.Color = 'red'; t_a.FontSize = 14;
    
    text_b = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_b_max, P_r_b_max/1e6, V_r_b_max/1e3);
    t_b = text(40, 20, text_b); t_b.Color = 'blue'; t_b.FontSize = 14;
    
    text_c = sprintf("max @ |Z_load| = %0.2f\n%0.2f MW, %0.2f kV",... 
        z_load_c_max, P_r_c_max/1e6, V_r_c_max/1e3);
    t_c = text(140, 50, text_c); t_c.Color = 'green'; t_c.FontSize = 14;
    %
    %
    fprintf("\n\n\t\tPowers in (MW), Voltages in (kV)\n");
    fprintf("|Z_load|  :     V_r_0  :     P_r_0   :     V_r_30  :     P_r_30  :    V_r_-30  :   P_r_-30\n");
    fprintf("-------------------------------------------------------------------------------------------\n");
    for ii = 1:1:size_z
        fprintf("%6d    :  %9.3f :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f\n",...
            z_load_absz(ii),...
            V_r_a(1, ii) / 1e3, P_r_a(ii) / 1e6,...
            V_r_b(1, ii) / 1e3, P_r_b(ii) / 1e6,...
            V_r_c(1, ii) / 1e3, P_r_c(ii) / 1e6...
            );
    end
    %
    %
    fprintf("\n\n\t\tPowers in (MW), Voltages in (kV), voltage reduced 5%%\n");
    fprintf("|Z_load|  :     V_r_0  :     P_r_0   :     V_r_30  :     P_r_30  :    V_r_-30  :   P_r_-30\n");
    fprintf("-------------------------------------------------------------------------------------------\n");
    for ii = 1:1:size_z
        fprintf("%6d    :  %9.3f :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f\n",...
            z_load_absz(ii),...
            V_r_a_down5(1, ii) / 1e3, P_r_a_down5(ii) / 1e6,...
            V_r_b_down5(1, ii) / 1e3, P_r_b_down5(ii) / 1e6,...
            V_r_c_down5(1, ii) / 1e3, P_r_c_down5(ii) / 1e6...
            );
    end
    %
    %
    fprintf("\n\n\t\tPowers in (MW), Voltages in (kV), voltage increased 5%%\n");
    fprintf("|Z_load|  :     V_r_0  :     P_r_0   :     V_r_30  :     P_r_30  :    V_r_-30  :   P_r_-30\n");
    fprintf("-------------------------------------------------------------------------------------------\n");
    for ii = 1:1:size_z
        fprintf("%6d    :  %9.3f :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f  :  %9.3f\n",...
            z_load_absz(ii),...
            V_r_a_up5(1, ii) / 1e3, P_r_a_up5(ii) / 1e6,...
            V_r_b_up5(1, ii) / 1e3, P_r_b_up5(ii) / 1e6,...
            V_r_c_up5(1, ii) / 1e3, P_r_c_up5(ii) / 1e6...
            );
    end
    %
end


%%%%%%%%~~~~~~~~~END>
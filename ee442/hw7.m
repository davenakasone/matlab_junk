%{
    hw7
    
    1  : fucking around
    2  : ch32, p11, buck
    3  : ch32, p12, boost
%}
clc;
clear;
close all;
select = 2;

if select == 1
    V_bgr = 1.25;    % bandgap reference voltage == V_out * R_B / (R_A + R_B)
    V_out = 2;       % must == V_bgr * (R_A + R_B) / R_B
    V_s = 5;         % source voltage
    D = 0.4;         % duty cycle, on-time portion
    f = 10e6;        % switching frequency

    R_load = 40;          % to simulate load on sps
    R_A = 6e3;            % divider resistor on negative terminal of opamp
    R_B = 10e3;           % divider resistor on negative terminal of opamp
    R_2 = 206;            % integrating resistor on opamp
    L = 100e-6;           % inductor value, H
    C_filter = 7.5e-9;    % filter capacitor at end of sps
    C_1 = 7.7e-9;         % opamp intergrator capacitor 1 of 2

    T = 1/f;                   % switching period
    I_out = V_out / R_load;    % nominal load current demanded from sps
    R_1 = f_para(R_A, R_B);    % thevinin resistance on negative input terminal of opamp
    C_2 = C_1 / 10;            % opamp integrator capacitor 2 of 2, smaller than C_1

    % exact TF
    syms s;
    Hs_loop = (1 + s * R_2 * C_1) / (s * (s * R_2 * C_1 * C_2 + C_1 + C_2));
    H_tf_loop = f_Fs_2_tf(Hs_loop);
    %f_rlocus(1/(1+H_tf_loop));
    %f_bode(1/(1+H_tf_loop));
    %f_pzplot(1/(1+H_tf_loop));
    %f_nyquist(1/(1+H_tf_loop));
    aa = R_B / (R_A + R_B);
    bb = (R_A + R_B) / (R_A * R_B);
    cc = (s * C_1) / (1 + s * C_1 * R_2);
    dd = s * C_2;
    ee = 1 / R_A;
    xx = aa * (bb + cc + dd) - ee;
    yy = cc + dd;
    %Hs = f_Fs_2_tf(xx/yy)
    %f_rlocus(Hs);
    %f_bode(Hs);
    %f_pz(Hs);
    %f_nyquist(Hs);

    beta = R_B / (R_A + R_B);
    AOLs = (1 + s * C_1 * R_2) / (s * (C_1 + C_2) + (s^2) * (C_1 * C_2 * R_2));
    ACLs = AOLs / (1 + beta * AOLs);
    tf_AOL = f_Fs_2_tf(AOLs);
    tf_ACL = f_Fs_2_tf(ACLs);
    
    f_bode([tf_ACL, tf_AOL, tf_ACL+tf_AOL]);
    f_pz([tf_ACL, tf_AOL]);
    f_rlocus([tf_ACL, tf_AOL]);
    f_nyquist([tf_ACL, tf_AOL, tf_ACL+tf_AOL]);
    f_step([tf_ACL, tf_AOL, tf_ACL+tf_AOL]);
    f_impulse([tf_ACL, tf_AOL, tf_ACL+tf_AOL]);
end


if select == 2
    repz = 3;
    cstep = [40e-9, 10e-9, 5e-9];
    
    V_bgr = 1.25;       % bandgap reference voltage == V_out * R_B / (R_A + R_B)
    V_out = 2;          % must == V_bgr * (R_A + R_B) / R_B
    V_s = 5;            % source voltage
    D = 0.4;            % duty cycle, on-time portion
    f_switch = 10e6;    % switching frequency

    R_load = 40;          % to simulate load on sps
    R_A = 6e3;            % divider resistor on negative terminal of opamp
    R_B = 10e3;           % divider resistor on negative terminal of opamp
    R_2 = 206;            % integrating resistor on opamp
    L = 100e-6;           % inductor value, H
    C_filter = 7.5e-9;    % filter capacitor at end of sps
    C1 = 7.7e-9;          % opamp intergrator capacitor 1 of 2
    C_1 = zeros(1, repz);        

    T_switch = 1/f_switch;     % switching period
    I_out = V_out / R_load;    % nominal load current demanded from sps
    R_1 = f_para(R_A, R_B);    % thevinin resistance on negative input terminal of opamp
    C2 = C1/10;                % opamp integrator capacitor 2 of 2, smaller than C_1
    C_2 = zeros(1, repz);    
    
    syms s;
    beta = R_B / (R_A + R_B);
    f_unity = zeros(1, repz);
    f_zero = zeros(1, repz);
    fprintf("\nV_bgr     = V_out * (R_B / (R_A + R_B) =  %0.3f V\n", V_out * (R_B/ (R_A + R_B)));
    fprintf("\n\tf_unity < f_zero/10 < f_resonant ?\n\n");
    for xx = 1:1:repz
        C_1(1,xx) = cstep(1,xx);
        C_2(1,xx) = C_1(1,xx)/10;
        
        fprintf("\nC1= %0.3f nF, C_2 = %0.3f nF\n", C_1(1,xx)*1e9, C_2(1,xx)*1e9);
        f_unity(1,xx) = 1 / (2 * pi * R_1  * C_1(1,xx));% unity or cross-over, == 0 dBm "|H| == 1"
        fprintf("f_uninty  =  %0.3f kHz\n", f_unity(1,xx)/1e3);
        f_zero(1,xx) = 1 / (2 * pi * R_2 * C_1(1,xx));% frequency of the zero, |H| == 0
        fprintf("f_zero    =  %0.3f kHz,   f_zero/10=  %0.3f  KHz\n",...
            f_zero(1,xx)/1e3, f_zero(1,xx)/1e4);
        f_resonant = 1 / (2 * pi * sqrt(L * C_filter));    % resonant frequency
        fprintf("f_resonant=  %0.3f  kHz\n", f_resonant/1e3);
        fprintf("f_switch  =  %0.3f MHz, T_switch=  %0.3f ns\n", f_switch/1e6, T_switch*1e9);
        
        ACLs = -1 * (1 + s * C_1(1,xx) * R_2) / (s * R_1 * C_1(1,xx));
        tfs_a(1,xx) = f_Fs_2_tf(ACLs);
        
        theo_AOLs = (1 + s * C_1(1,xx) * R_2) /...
            (s * (C_1(1,xx) + C_2(1,xx)) + (s^2) * (C_1(1,xx) * C_2(1,xx) * R_2));
        tfs_b(1,xx) = f_Fs_2_tf(beta * theo_AOLs);
        theo_ACLs = 1 / (1 + beta * theo_AOLs);
        tfs_c(1,xx) = f_Fs_2_tf(theo_ACLs);

        aa = R_B / (R_A + R_B);
        bb = (R_A + R_B) / (R_A * R_B);
        cc = s * C_1(1,xx) / (1 + s * R_2 * C_1(1,xx));
        dd = s * C_2(1,xx);
        ee = 1 / R_A;
        ff = aa * (bb + cc + dd) - ee;
        gg = cc + dd;
        GHs = ff/gg;
        tfs_d(1,xx) = f_Fs_2_tf(GHs);
    end
    f_bode([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_pz([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);
    %{
    f_bode([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_bode([tfs_b(1,1), tfs_b(1,2), tfs_b(1,3)]);
    f_bode([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);

    f_rlocus([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_nyquist([tfs_b(1,1), tfs_b(1,2), tfs_b(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);    
    f_impulse([tfs_d(1,1), tfs_d(1,2), tfs_d(1,3)]);  
    %}
end


if select == 3
    repz = 3;
    cstep = [300e-9, 50e-9, 10e-9];
    
    V_bgr = 1.25;       % bandgap reference voltage == V_out * R_B / (R_A + R_B)
    V_out = 15;          % must == V_bgr * (R_A + R_B) / R_B
    V_s = 5;            % source voltage
    D = 2/3;            % duty cycle, on-time portion
    f_switch = 1e6;    % switching frequency

    R_load = 40;          % to simulate load on sps
    R_A = 110e3;            % divider resistor on negative terminal of opamp
    R_B = 10e3;           % divider resistor on negative terminal of opamp
    R_2 = 183;            % integrating resistor on opamp
    L = 4700e-6;           % inductor value, H
    C_filter = 0.22e-6;    % filter capacitor at end of sps
    C1 = 0.174e-6;          % opamp intergrator capacitor 1 of 2
    C_1 = zeros(1, repz);        

    T_switch = 1/f_switch;     % switching period
    I_out = V_out / R_load;    % nominal load current demanded from sps
    R_1 = f_para(R_A, R_B);    % thevinin resistance on negative input terminal of opamp
    C2 = C1/10;                % opamp integrator capacitor 2 of 2, smaller than C_1
    C_2 = zeros(1, repz);    
    
    syms s;
    beta = R_B / (R_A + R_B);
    f_unity = zeros(1, repz);
    f_zero = zeros(1, repz);
    fprintf("\nV_bgr     = V_out * (R_B / (R_A + R_B) =  %0.3f V\n", V_out * (R_B/ (R_A + R_B)));
    fprintf("\n\tf_unity < f_zero/10 < f_resonant ?\n\n");
    for xx = 1:1:repz
        C_1(1,xx) = cstep(1,xx);
        C_2(1,xx) = C_1(1,xx)/10;
        
        fprintf("\nC1= %0.3f nF, C_2 = %0.3f nF\n", C_1(1,xx)*1e9, C_2(1,xx)*1e9);
        f_unity(1,xx) = 1 / (2 * pi * R_1  * C_1(1,xx));% unity or cross-over, == 0 dBm "|H| == 1"
        fprintf("f_uninty  =  %0.3f kHz\n", f_unity(1,xx)/1e3);
        f_zero(1,xx) = 1 / (2 * pi * R_2 * C_1(1,xx));% frequency of the zero, |H| == 0
        fprintf("f_zero    =  %0.3f kHz,   f_zero/10=  %0.3f  KHz\n",...
            f_zero(1,xx)/1e3, f_zero(1,xx)/1e4);
        f_resonant = 1 / (2 * pi * sqrt(L * C_filter));    % resonant frequency
        fprintf("f_resonant=  %0.3f  kHz\n", f_resonant/1e3);
        fprintf("f_switch  =  %0.3f MHz, T_switch=  %0.3f ns\n", f_switch/1e6, T_switch*1e9);
        
        ACLs = -1 * (1 + s * C_1(1,xx) * R_2) / (s * R_1 * C_1(1,xx));
        tfs_a(1,xx) = f_Fs_2_tf(ACLs);
        
        theo_AOLs = (1 + s * C_1(1,xx) * R_2) /...
            (s * (C_1(1,xx) + C_2(1,xx)) + (s^2) * (C_1(1,xx) * C_2(1,xx) * R_2));
        tfs_b(1,xx) = f_Fs_2_tf(beta * theo_AOLs);
        theo_ACLs = 1 / (1 + beta * theo_AOLs);
        tfs_c(1,xx) = f_Fs_2_tf(theo_ACLs);

        aa = R_B / (R_A + R_B);
        bb = (R_A + R_B) / (R_A * R_B);
        cc = s * C_1(1,xx) / (1 + s * R_2 * C_1(1,xx));
        dd = s * C_2(1,xx);
        ee = 1 / R_A;
        ff = aa * (bb + cc + dd) - ee;
        gg = cc + dd;
        GHs = ff/gg;
        tfs_d(1,xx) = f_Fs_2_tf(GHs);
    end
    f_bode([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_pz([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);
    %{
    f_bode([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_bode([tfs_b(1,1), tfs_b(1,2), tfs_b(1,3)]);
    f_bode([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);

    f_rlocus([tfs_a(1,1), tfs_a(1,2), tfs_a(1,3)]);
    f_nyquist([tfs_b(1,1), tfs_b(1,2), tfs_b(1,3)]);
    f_step([tfs_c(1,1), tfs_c(1,2), tfs_c(1,3)]);    
    f_impulse([tfs_d(1,1), tfs_d(1,2), tfs_d(1,3)]);  
    %}
end


if select == 99
    fprintf("\ndone\n");
end


%%%%%%%%~~~~~~~~END>  hw7.m

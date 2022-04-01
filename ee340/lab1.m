%{
    lab1

    1  : question 1a, impedance and phase angle
    2  : question 1b, find equiv parallel
    3  : question 1c, add capacitor
    4  : question 1d, put capacitor in for unity (all real power)
%}
clc;
select = 4;

%-------------------------------------------------------------------------------------
if select == 1
    V_s = 120;                             % volts supplied
    f = 60;                                % Hz 
    w = 2 * pi * f;                        % omega as rad/sec
    w_draw_motor = 192;                    % watts used by motor
    w_produce_motor = f_hp2watts(0.25);    % watts produced, convert 1/4 hp {186.42 W}
    i_motor = 2;                           % amps used by motor
    
    S_mag = abs(V_s) * abs(i_motor);
    fprintf("apparent power |S| = |V||I| = %d * %d = %d VA\n",...
        V_s, i_motor, S_mag);
    
    P = w_draw_motor;
    Q = sqrt(S_mag^2 - P^2);
    fprintf("\nQ = %d VAR\n", Q);
    
    pf = P / S_mag;
    fprintf("\npf = P / |S| = %0.2f\n", pf);
    tht = acos(pf);
    fprintf("\nangle = %0.3f  rad  ,  %0.3f  deg\n", tht, rad2deg(tht));
    
    Z_mag = abs(V_s) / abs(i_motor);
    R = w_draw_motor / i_motor^2;
    X = sqrt(Z_mag^2 - R^2);
    Z = R + 1j*X;
    fprintf("\n");
    f_mdri("Z", Z, 1);
end

%-------------------------------------------------------------------------------------
if select == 2
    V_s = 120;                             % volts supplied
    f = 60;                                % Hz 
    w = 2 * pi * f;                        % omega as rad/sec
    w_draw_motor = 192;                    % watts used by motor
    w_produce_motor = f_hp2watts(0.25);    % watts produced, convert 1/4 hp {186.42 W}
    i_motor = 2;                           % amps used by motor
    
    S_mag = abs(V_s) * abs(i_motor);
    P = w_draw_motor;
    Q = sqrt(S_mag^2 - P^2);
    pf = P / S_mag;
    tht = acos(pf);
    
    Z_mag = abs(V_s) / abs(i_motor);
    R = w_draw_motor / i_motor^2;
    X = sqrt(Z_mag^2 - R^2);
    Z = R + 1j*X;
    
    R_new = V_s^2 / P;                  % 75
    X_new = V_s^2 / Q;                  % 100
    check = f_para(R_new, 1j*X_new);    % 48.0000 +36.0000i ...same Z
end

%-------------------------------------------------------------------------------------
if select == 3
    V_s = 120;                             % volts supplied
    f = 60;                                % Hz 
    w = 2 * pi * f;                        % omega as rad/sec
    w_draw_motor = 192;                    % watts used by motor
    w_produce_motor = f_hp2watts(0.25);    % watts produced, convert 1/4 hp {186.42 W}
    i_motor = 2;                           % amps used by motor
    
    S_mag = abs(V_s) * abs(i_motor);
    P = w_draw_motor;
    Q = sqrt(S_mag^2 - P^2);
    pf = P / S_mag;
    tht = acos(pf);
    
    X_c = 200;                          % capacitor in parallel
    Q_c = V_s^2 / X_c;                  % 72 var, lowers total Q
    Q_new = Q - Q_c;                    % 144 - 72 = 72
    S_mag_new = sqrt(P^2 + Q_new^2);    % 205.1 VA
    i_new = S_mag_new / abs(V_s);       % 1.71 A
    pf_new = P / S_mag_new              % pf = 0.936
    
    S_new = P + 1j*Q_new;
    check = angle(S_new);
    checkk = cos(check);
end

%-------------------------------------------------------------------------------------
if select == 4
    V_s = 120;                             % volts supplied
    f = 60;                                % Hz 
    w = 2 * pi * f;                        % omega as rad/sec
    w_draw_motor = 192;                    % watts used by motor
    w_produce_motor = f_hp2watts(0.25);    % watts produced, convert 1/4 hp {186.42 W}
    i_motor = 2;                           % amps used by motor
    
    S_mag = abs(V_s) * abs(i_motor);
    P = w_draw_motor;
    Q = sqrt(S_mag^2 - P^2);
    
    X_c = V_s^2 / Q;                % 100, need it to cancel inductor
    S_mag_new = P;                  % 192, because P = |S| pf = |S| 1
    Q_new = 0;                      % because unity
    i_new = S_mag_new / abs(V_s);   % 1.6 A
end

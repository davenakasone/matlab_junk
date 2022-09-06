%{
    lab final 

    1  :  power Q of an inductor
    2  :  capacitor bank
    3  :  power factor
    4  :  generator
    5  :  3-phase induction motor 
    6  :  an electric bill
    7  :  
    
%}
close all;
clear;
clc;
select = 7;


%-------------------------------------------------------------------------------------
if select == 1
    % inductor on an AC source
    f = 60;
    w = 2 * pi * f;
    V = 120;
    L = 50e-3;
    z = 1j * w * L;
    I = V / z;
    S = V * conj(I);
    fprintf("it draws:  %0.1f VAR\n", imag(S));
end


%-------------------------------------------------------------------------------------
if select == 2
    V_rate = 240;
    f_rate = 60;
    w_rate = 2 * pi * f_rate;
    Q_rate = 30e3;
    I_rate = (Q_rate / V_rate);
    X_rate = Q_rate / abs(I_rate)^2;
    C = 1 / (X_rate * w_rate);

    f = 50;
    w = 2 * pi * f;
    V = 230;
    z = 1 / (1j * C * w);
    I = V / z;
    S = V * conj(I);
    fprintf("Q=  %0.2f kVAR   ...supplied\n", abs(imag(S))/1e3);
end


%-------------------------------------------------------------------------------------
if select == 3
    R = 4;
    X = 2;
    z = R + 1j*X;   % must lag...
    theta = angle(z);
    pf = cos(theta);
    fprintf("pf =  %0.1f %%  ,  angle=  %0.1f deg\n", 100*pf, rad2deg(theta));
end


%-------------------------------------------------------------------------------------
if select == 4
    V_rate = 25e3;
    S_rate = 50e6;
    z_base = V_rate^2 / S_rate;
    Xpu = 0.05;
    X = z_base * Xpu;
    fprintf("it is actually:  %0.3f  ohm\n", X);
end


%-------------------------------------------------------------------------------------
if select == 5
    hp_rate = 10;
    P_hp = f_hp2watts(hp_rate);
    V_rate = 480;
    I_rate = 10;
    pf = 0.95;
    theta = acos(pf);
    
    P_in = sqrt(3) * V_rate * I_rate * pf;
    P_out = P_hp;
    eff_eta = P_out / P_in;
    fprintf("eff =  %0.1f %%\n", 100*eff_eta);
end


%-------------------------------------------------------------------------------------
if select == 6
    bill = 135.72;
    rate = 0.13; % per kWh
    consume = bill / (rate * 24 * 30); % hours in day * days in month
    fprintf("avg consumption:  %0.2f kW\n", consume);
end


%-------------------------------------------------------------------------------------
if select == 7
    
end


%-------------------------------------------------------------------------------------
if select == 99
 
end


%%%%%%%%~~~~~~~~~END> 

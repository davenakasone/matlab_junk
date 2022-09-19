%{
    v2, ch2, kenetic theory of gas
    
    0  :  testing, ex 2.10
    1  :  q1, gas concepts
    2  :  q2, find moles
    3  :  q3, ideal gas, concepts 
    4  :  q4, gauge pressure
    5  :  q5, bike tire
    6  :  q6, gas in cylinder
    7  :  q7, particle density
    8  :  q8, hot air ballon
    9  :  q9, speeds
    10 :  q10, kentic energy
    11 :  q11, escape velocity
    12 :  q12, fusion temperatures
    13 :  q13, vrms and temperature
    14 :  q14, uranium gas diffusion
    15 :  q15, work

%}
clc;
clear;
select = 15;


if select == 0
    M_nitrogren = 28e-3; % kg/mol
    T = konst.celsius_2_kelvin(27); % K
    rat = konst.max_bolt_v(300, 4.65e-26, T) / konst.max_bolt_v(100, 4.65e-26, T)
end


%%%%%~~~~


if select == 1
    fprintf("a) P' = (1/6)P\n");
    fprintf("b) P'' = (1/3)P\n");
end


%%%%%~~~~


if select == 2
    V = 4.7 / 1000; % m^3 of air
    P = konst.atm_2_pa * 1; % Pa
    T = 35.5 + 273.15;
    n = P * V / (8.314 * T); % moles
    fprintf("moles:  %0.4f\n", n);
end


%%%%%~~~~


if select == 3
    fprintf("a)  T = P / (rho k_B)\n");

    rho = 0.95e27;
    P = 89;
    T = P * (1.01e5) / (rho * 1.38e-23);
    fprintf("b)  T =  %0.3f C\n", T-273.15);

    r = 5.7e15;  % m
    m = 4000 * (1.989e30);    % times sun
    Pn = 6.9e-9; % Pa
    temp = konst.volume_of_sphere(r) * 2 * Pn * (1.67e-27);
    T = temp / (m * (1.38e-23));
    fprintf("c)  T =  %0.3f K\n", T);

    rho_0 = 1.08 * 10^6; % particle density, space
    P0 = 4.8e-17; % N/m^2
    T0 = P0 / (rho_0 * (1.38e-23));
    fprintf("4)  T0 =  %0.3f K\n", T0);

end


%%%%%~~~~


if select == 4
    P1 = 3.5e5; % Pa
    T1 = 35 + 273.15; % K
    T2 = -40 + 273.15; % K
    P2 = P1 * (T2/T1)
    fprintf("a)  P =  %0.3f  Pa\n", P2);
end


%%%%%~~~~


if select == 5
    P1 = 6.9e5; % Pa
    T = 19.5 + 273.15;
    V = 2 / 1000;
    Vout = 85 / 10^6;

    N1 = P1 * V / ((1.38e-23) * T);           % molecules before
    dN = (1.013e5) * Vout / ((1.38e-23) * T); % molecules out
    N2 = N1 - dN;
    P2 = N2 * (1.38e-23) * T / V;
    fprintf("P2 = %0.2f x 10^5  Pa\n", P2/10^5);
end


%%%%%~~~~


if select == 6
    Vg = 50 / 1000;
    P1 = 1.45e7;
    T1 = 25 + 273.15;
    T2 = -78.5 + 273.15;

    P2 = P1 * (T2/T1);
    fprintf("a)  P2 =  %0.2f x 10^6 Pa\n", P2/10^6);

    P2 = P1 * (T2/T1) * 0.9;
    fprintf("b)  P2 =  %0.2f x 10^6 Pa\n", P2/10^6);

    T = T1 * ((1.013e5)/P1);
    fprintf("c)  T = %0.2f K\n", T);
    fprintf("d) fuck no\n");

end


%%%%%~~~~


if select == 7
    rhoN = 10^6; % atoms/m^3
    T = 2.7; % K

    P = rhoN * (1.38e-23) * T;
    fprintf("a) pressure = %0.1f e-17 Pa\n", P*10^17);

    V = 2.5 * 8.31 * T / P;
    fprintf("b) V =  %0.1f e17 m^3\n", V/10^17);

    L = V^(1/3) / 1000;
    fprintf("c)  L = %0.1f km\n", L);
end


%%%%%~~~~


if select == 8
    Vmax = 1920; % m^3
    mp = 390; % kg
    M_air = 28.97 / 1000; % kg/mol
    rho_air = 1.2; % kg/m^3
    
    Ta = 101325 * Vmax * M_air;
    Tb = 8.31 * (rho_air * Vmax - mp);
    T = (Ta/Tb) - 273.15;
    fprintf("a)  T = %0.3f C\n", T);

    Tc = 101325 * Vmax * M_air * 0.8;
    Td = 8.31 * (rho_air * Vmax * 0.85 - mp) ;
    Tmax = (Tc/Td) - 273.15;
    fprintf("b)  T = %0.3f C\n", Tmax);
end


%%%%%~~~~


if select == 9
    m = 6.63e-26; % kg
    T = 2550; % K
    vrms = sqrt((3 * (1.38e-23) * T) / m);
    fprintf("a) vrms:  %0.3f m/s\n", vrms);


end


%%%%%~~~~


if select == 10
    T1 = 11 + 273.15; % K
    T2 = 22 + 273.15; % K
    
    Trat = T2/T1;
    fprintf("ratio:  %0.3f  infer percent   about the same\n", Trat);
end


%%%%%~~~~


if select == 11
    M = 2.016 / 1000; % kg/mol
    vesp = 2.38 * 1000; % m/s
    
    m = M / (6.02e23); % kg
    T = m * vesp^2 / (3 * (1.38e-23));
    fprintf("T= %0.1f\n", T);
end


%%%%%~~~~


if select == 12
    E = 5.8e-14; % J/ atom
    T = 2 * E / (3 * (1.38e-23));
    fprintf("T= %0.2f e9 K\n", T/10^9);
end


%%%%%~~~~


if select == 13
    M = 44/1000; % kg/mol
    m = M / (6.02e23);
    vrms = 1.1e5;
    T = m * vrms^2 / (3 * (1.38e-23));
    fprintf("T= %0.2f e7 K\n", T/10^7);
end

%%%%%~~~~


if select == 14
    M235 = 349/1000; % kg/mol
    M238 = 352/1000; % kg/mol

    Urat = sqrt(M238/M235);
    fprintf("a) ratio:  %0.3f\n", Urat);
    fprintf("\nv235 = %0.3f * v238 = v238 + %0.3f * v238\n", Urat, Urat-1);
    fprintf("dv = %.03f * v238\n", Urat-1);

    m238 = M238/(6.02e23);
    dv = Urat-1;
    temp = 3 * (1.38e-23);
    T = m238 * (1/dv)^2 / temp;
    fprintf("c) T = %0.1f K   ...%0.1f F\n", T, konst.kelvin_2_farhrenheit(T));
end


%%%%%~~~~


if select == 15
    work = 1100; % J
    Pf = 1e5; % Pa
    V = 0.027; % m^3
    R = 8.314;

    Vf = V * exp(-1 * work / (Pf * V));
    fprintf("Vf = %0.5f m^3\n", Vf);

    T = Pf * V / R;
    fprintf("T = %0.3f K\n", T);

end


%%%%%~~~~


if select == 99
    fprintf("\n\tDONE\n");
end


%%%%%%%%~~~~~~~~END>  v1_ch14_flud.m
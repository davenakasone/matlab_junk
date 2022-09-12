%{
    volume 2, chapter 1, heat and temperature

    0   :  try conversions
    -1  :  some book work
    -2  :  ex1.13
    1   :  p1, temperature
    2   :  p2, thermal expansion, use table alpha
    3   :  p3, thermal expansion, use table beta, then divide
    4   :  p4, overflow from expansion, using beta
    5   :  p5, density, pressure, states
    6   :  p6, water states, some chemistry, conserve energy, bullshit problem
    7   :  p7, specific heat
    8   :  p8, conserve energy
    9   :  p9, heat transfer + phase change
    10  :  p10, iceberg melts
    11  :  p11, window conducts
    12  :  p12, radiation


%}
clc;
clear;
select = 12;


if select == 0
    temp_a = 0;
    temp_b = konst.celsius_2_fahrenheit(temp_a);
    fprintf("c2f  :  %0.3f C -->  %0.3f F\n", temp_a, temp_b);
    temp_c = konst.celsius_2_kelvin(temp_a);
    fprintf("c2k  :  %0.3f C -->  %0.3f K\n", temp_a, temp_c);
    temp_d = konst.kelvin_2_celsius(temp_c);
    fprintf("k2c  :  %0.3f K  -->  %0.3f C\n", temp_c, temp_d);
    temp_g = konst.kelvin_2_farhrenheit(temp_c);
    fprintf("k2f  :  %0.3f K  -->  %0.3f F\n", temp_c, temp_g)
    temp_e = konst.fahrenheit_2_celsius(temp_b);
    fprintf("f2c  :  %0.3f F  -->  %0.3f C\n", temp_b, temp_e);
    temp_f = konst.fahrenheit_2_kelvin(temp_b);
    fprintf("f2k  :  %0.3f F  -->  %0.3f K\n", temp_b, temp_f);
end


%%%%%~~~~


if select == -1
    % 2 rods, aluminum, steel, same size, put isoloated in diff temp, welded end to end
    r = (1/100) / 2; % m
    A = konst.area_of_circle(r); % m^2
    L = 0.25; % m
    k_alu = 220; % W/(m C)
    k_stl = 80; % W/(m C)
    T_stl = 100; % C
    T_alu = 20; % C

    syms Tf;
    Pstl = (k_stl * A * (T_stl - Tf)) / L;
    Palu = -1*(k_alu * A * (T_alu - Tf)) / L;
    T_f = vpasolve(Pstl==Palu, Tf);
    P_stl = subs(Pstl, Tf, T_f);
    P_alu = subs(Palu, Tf, T_f);
    fprintf("power for each:  %0.3f  %0.3f   W\n", P_stl, P_alu);
    
    d_T = T_stl - T_alu; % C
    % see p52, ch1, v1
    % can do it differently with R values... R=L/k 
end


%%%%%~~~~


if select == -2
    e = 0.97; % emisstivity
    T_obj = konst.celsius_2_kelvin(33);
    T_env = konst.celsius_2_kelvin(22);
    A = 1.5; % m^2
    P_net = konst.boltzmann_radiation * e * A * (T_env^4 - T_obj^4);
    fprintf("%0.3f\n", P_net);
end


%%%%%~~~~


if select == 1
    t_c = 49; % C

    fprintf("a) it should be in a reasonable range...\n");
    fprintf("b) c2f  :  %0.1f F\n", konst.celsius_2_fahrenheit(t_c));
    fprintf("c) it should be in a reasoable range...\n");
    fprintf("d) c2k  :  %0.1f K  no decimal\n", konst.celsius_2_kelvin(t_c));
end


%%%%%~~~~


if select == 2
    % height measured, at temperature, then changes...
    L_0 = 170; % m
    T_0 = 35; % C
    T_1 = -10; % C

    d_T = T_1 - T_0; % C
    alpha_marble = 2.5e-6; % see table, thermal coeff
    d_L = alpha_marble * d_T;
    L = L_0 + d_L;
    fprintf("the mountain is now:  %0.5f m    169.98 is what they want\n", L);
end


%%%%%~~~~


if select == 3
    % global warming, ocean uses beta, but /3 for rise...
    % change in length for column of water 1km high, 1 C
    L_0 = 1e3; % m
    d_C = 1; % C
    beta_water = 210e-6;
    
    d_L = (beta_water/3) * L_0 * d_C;
    fprintf("change in level is:  %0.3f m   ...use 7e-3? yes\n", d_L);
end


%%%%%~~~~


if select == 4
    % radiatior may overflow when hot
    V_l = 14; % L  radiator resivor made of copper
    T_0 = 10; % C
    T_1 = 95; % C, operating temperature
    beta_cool = 400e-6;
    beta_cop = 5.1e-5;

    dV = (beta_cool - beta_cop) * V_l * (T_1 - T_0);
    fprintf("overflow is:  %0.3f L\n", dV);
end


%%%%%~~~~


if select == 5
    rho_water_liquid = 1000; % kg/m^3
    rho_water_ice = 917; % kg/m^3

    dV = rho_water_liquid / rho_water_ice; % expansion relative to ice
    pressure = ((dV-1)/1) * (2.2e9);
    fprintf("a) pressure=  %0.2f  x 10^8  N/m^2     use 1.98e8\n", pressure/1e8);
    fprintf("b) big pressure, give cells room to expand\n");
end


%%%%%~~~~


if select == 6
    % water evaporated should == energy when it condenses...
    % mole of water is 2H + 1O
    rain_volume = (4.9e5)*(1e9); % m^3   
    avg_h = 8.4 * (1e3); % m
    J_evap_mole = 40.6e3; % J  of heat per mole
    rho_water = 1000;

    m_water_amu = 2 * 1.0079 + 15.9994;
    m_water_amu_kg = m_water_amu * (1.6605e-27);
    m_water_kg_per_mol = m_water_amu_kg * (6.02214179e23);
    E_kg = J_evap_mole / m_water_kg_per_mol;
    mass = rho_water * rain_volume;

    E_evap = mass * E_kg;
    E_cond = mass * 9.8 * avg_h;
    Eng = (E_evap + E_cond) / 365.25
end


%%%%%~~~~


if select == 7
    T_env = 0; % C
    breath_per_min = 17;
    inhale_V = 0.00045; % m^3 of air
    T_obj = 37; % C
    rho_air = 1.2; % kg/m^3
    c_air = 1e3; % J/(kg C)
    
    fprintf("a) energy to warm 1 breath:  Q = mc dT = mc(Tbody - Tenv)\n");
    mass_air = inhale_V * rho_air;
    dT = T_obj - T_env;
    Q = mass_air * c_air * dT;
    breath_per_sec = breath_per_min / 60;
    P = Q * breath_per_sec;
    fprintf("b)  P=  %0.3f W\n", P);
end


%%%%%~~~~


if select == 8
    % wood moves, converts all to energy
    c_wood = 130; % J/(kg C)
    h = 1.5; % m
    fprintf("a) see equations\n");

    dT = (9.81 * h) / c_wood;
    fprintf("K==C,  dT=  %0.5f\n", dT);
end


%%%%%~~~~


if select == 9
    m_ice = 0.22; % kg
    T_i = -20; % C
    T_f = 130; % C
    
    % heat it 0 C
    dT_0 = 0 - T_i;
    c_ice = 2.09; % kJ/(kg C)
    Q_0 = round(m_ice * c_ice * dT_0, 3, "significant") % kJ

    % melt it at 0 C
    Lf = 334; % kJ/kg
    Q_1 = round(m_ice * Lf, 3, "significant") % kJ

    % warm it to 100 C
    dT_1 = 100 - 0;
    c_water = 4.186; % kJ/(kg C)
    Q_2 = round(m_ice * c_water * dT_1, 3, "significant")

    % vaporize it at 100 C
    Lv = 2256; % kJ/kg
    Q_3 = round(m_ice * Lv, 3, "significant")

    % heat the vapor...
    dT_2 = T_f - 100;
    c_vapor = 1.52; % kj/(kg C)
    Q_4 = round(m_ice * c_vapor * dT_2, 3, "significant");

    Q_kj = Q_0 + Q_1 + Q_2 + Q_3 + Q_4;
    Q_kcal = round(konst.joule_2_kcal * 1000 * Q_kj, 3, "significant");
    fprintf("a)  E =  %0.1f kcal\n", Q_kcal);

    % time at each stage,  20.0 kJ/s
    P = 20.0; % kJ/s
    t0 = round(Q_0 / P, 3, "significant")
    t1 = round(Q_1 / P, 3, "significant")
    t2 = round(Q_2 / P, 3, "significant")
    t3 = round(Q_3 / P, 3, "significant")
    t4 = round(Q_4 / P, 3, "significant")
    t = Q_kj / P

end


%%%%%~~~~


if select == 10
    L = 160e3;
    W = 40e3;
    H = 250;
    V = L*W*H;
    rho_ice = 917; % kg/m^3

    %mass = round(rho_ice * V, 3, "significant") % kg
    mass = rho_ice * V

    % heat in J to melt...
    Lf = 79.8 * 4186; % kcal -> J/kg
    Q_melt = round(mass * Lf, 3, "significant")  % J
    A = L*W;
    P = 100; % W/m^2
    %Q_day = round(P * A * 12 * 3600, 4, "significant") % J
    Q_day = P * A * 12 * 3600
    n = (Q_melt / Q_day) * (1/365.25)
end


%%%%%~~~~


if select == 11
    A = 1.25; % m^2
    d_glass = 0.745/100; % m, 2 panes
    d_air = 0.85/100; % m, air gap
    T_in = 15; % C
    T_out = -10; % C
    k_air = 0.023;
    k_glass = 0.84;

    dT = T_in - T_out;
    P = (k_glass * k_air * A * dT) / (2 * k_air * d_glass + k_glass * d_air)

    d = 1.46/100;
    singp = k_glass * A * dT / d
end


%%%%%~~~~


if select == 12
    T_obj = 101.5; % C
    T_env = 41.5; % C
    em = 0.75;
    A = 1.05; % m^2

    P = A * em * (5.67e-8) * ((T_obj+273.15)^4 - (T_env + 273.15)^4)
end


%%%%%~~~~


if select == 99
    fprintf("\n\tDONE\n");
end


%%%%%%%%~~~~~~~~END>  v2_ch1_heat_temp.m
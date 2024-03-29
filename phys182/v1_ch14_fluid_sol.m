%{
    volume 1, chapter 14, fluids

    1   :  hw1, pressure
    2   :  hw1, pressure
    3   :  hw1, pressure at depth
    4   :  hw1, pressure at depth
    5   :  hw1, pressure U tube
    6   :  hw1, pressure with density, ballon + bouyant
    7   :  hw1, pressure, bouyancy, floating
    8   :  hw1, archemides...
    9   :  hw1, pressure drop
    10  :  hw1, concrete through a hose
    11  :  hw1, flows
    12  :  hw1, flows with obscure units
    13  :  hw1, hoover dam
    14  :  hw1, airplane wing and pressure
    15  :  hw1, pressue in tank...

%}
clc;
clear;
select = 15;


if select == 1
    % weight momentarily placed on heel
    % find pressure, pounds per in^2
    area_cm2 = 1.25;                           % cm^2
    mass_kg = 51.5;                            % kg
    Fg_newton = mass_kg * konst.g;             % N,   kg (m/s^2)
    p_Pa_cm = Fg_newton / area_cm2;            % N/cm^2
    p_psi = p_Pa_cm * konst.newton_2_pound * konst.cm_per_inch^2;
    fprintf("the pressure is:  %0.3f psi\n", p_psi);
end


%%%%%~~~~


if select == 2
    % find out how much pressure teeth make
    area_mm2 = 0.75;                         % mm^2
    area_m2 = area_mm2 * (1/1000)^2;         % m^2
    force_newton = 475;                      % N, kg/m^2
    pressure_pa = force_newton / area_m2;    % Pa
    fprintf("the pressure is:  %0.3f Pa\n", pressure_pa);
end


%%%%%~~~~


if select == 3
    % basketball in a pool...doesn't change anything as it moves
    p_gauge = 55e3; % Pa, the gaugue pressure ball is pumped to
    p_atm = 101e3; % Pa, atmospheric pressure at the surface of the pool

    p_absolute = p_atm + p_gauge;   % pressure in ball at surface
    fprintf("a)  p_abs:  %0.1f kPa   ...no decimal\n", p_absolute/1000);

    fprintf("b) ball pressure stays the same as it sinks...no changes, so inside pressure SAME\n");
    fprintf("c) fluid pressure increases with depth...so outside pressure INCREASES\n");
    fprintf("d) [Delta]P = p_g - [rho] g y  ,  pressure difference as submereged distance y\n")

    y_p0 = p_gauge / (9.81 * konst.rho_water);
    fprintf("e)  no pressure difference at:  %0.3f m\n", y_p0);

    % ball is submerged in mercury now
    rho_mercurry = 13.5e3;
    y_p0_m = p_gauge / (9.81 * rho_mercurry);
    fprintf("f)  no pressure difference at:  %0.4f m\n", y_p0_m);
end


%%%%%~~~~


if select == 4
    % pressure at bottom of ocean
    rho_swater = 1.025e3; % kg/m^3
    deep = 11e3; % m
    p_pa = rho_swater * deep * konst.g;
    p_atmos = p_pa * konst.pa_2_atm;
    fprintf("pressure at bottom is:  %0.2f atm\n", p_atmos);
end


%%%%%~~~~


if select == 5
    % vertical U tube, 2 liquids don't mix
    fprintf("pressure on both sides must be equal\n");
    fprintf("p = [rho] g h    Pa\n");
    fprintf("p_left = p_right  ,  [rho_a] [g] [d1] = [rho_a] [g] [d2] + [rho_b] [g] [d3]\n");
    fprintf("solving for [rho_b],  [rho_b] = [rho_a] [d1-d2] / [d3]    kg/m^3\n");

    d1 = 10.1 / 100; % m
    d2 = 7.1 / 100; % m
    d3 = 6.1 / 100; % m
    rho_a = 1.5e3; % kg/m^3
    rho_b = rho_a * (d1-d2) / d3;
    fprintf("b)  [rho_b] =  %0.3f  kg/m^3\n", rho_b);
end


%%%%%~~~~


if select == 6
    % ballon is filled with helium, accelerates up
    rho_air = 1.23; % kg/m^3
    ballon_rho = 0.22; % kg/m^3
    ballon_volume = 0.021; % m^3

    % express height in terms of rho_he and volume
    ballon_mass = ballon_rho * ballon_volume;  % first find the mass
    ballon_weight = ballon_mass * 9.81; % now weight
    fprintf("a) ballon weight = density * volume * g\n");
    fprintf("b) ballon weight =  %0.5f  N\n", ballon_weight);

    % express bouyant force in terms of volume and density
    % remember it is how much of the other fluid the ballon displaces
    fprintf("c) F_bouy = [rho_air] * ballon_volume * g\n");
    f_bouy = rho_air * ballon_volume * 9.81;
    fprintf("d) F_bouy =  %0.4f N\n", f_bouy);

    % assume upward is positive, express total force
    fprintf("e) F_total = F_bouy - Weight    ...bouyant less gravity\n");
    f_total = f_bouy - ballon_weight;
    fprintf("f) F_total=  %0.4f N\n", f_total);

    % expressing acceleration in terms of F, ballon_rho, ballon_volume
    fprintf("g) a = F_total / ([rho_ballon] * [volume])   m/s^2\n");
    %a = f_total / (ballon_rho * ballon_volume); %don't use ?
    a = ((rho_air/ballon_rho) - 1) * 9.81;
    fprintf("h) a=  %0.3f  m/s^2\n", a);
end


%%%%%~~~~


if select == 7
    ball_radius = 11 / 100; % m
    ball_volume = (4/3)*pi*ball_radius^3;
    fprintf("3 balls have same radius: %0.3f m\n", ball_radius);
    fprintf("ball#1 floats, half exposed\n");
    b2_rho = 801; % kg/m^3
    fprintf("ball#2 anchored at bottom, has denisty: %0.1f kg/m^s\n", b2_rho);
    b3_rho = 1200; % kg/m^3
    fprintf("ball#3 held so fully submered...density:  %0.1f kg/m^3\n", b3_rho);
    rho_water = 1000; % kg/m^3
    fprintf("assume density of water is:  %0.1f  kg/m^3\n", rho_water);
    fprintf("a) floating ball is at rest, f_net == 0, abs(f_bouy) = weight\n");

    % tension on rope holding second ball
    % it would float...density is < water
    % sum(F) = F_bouy - mg - T = 0    it is in equalib
    b2_f_bouy = rho_water * ball_volume * 9.81;
    b2_weight = b2_rho * ball_volume * 9.81;
    b2_tension = b2_f_bouy - b2_weight;
    fprintf("b) ball#2 tension:  %0.3f N\n", b2_tension);

    % for 3rd ball, it is more dense..but held, find tension
    b3_f_bouy = rho_water * ball_volume * 9.81;
    b3_weight = b3_rho * ball_volume * 9.81;
    %b3_tension = b3_weight - b3_f_bouy; % reversed direction...
    b3_tension = (b3_rho-rho_water) * ball_volume * 9.81;
    fprintf("c) ball#3 tension:  %0.3f N\n", b3_tension);
end


%%%%%~~~~


if select == 8
    % can find density of fluid or a solid...    note g/cm^3
    iron_mass_air = 375; % g, mass in air
    iron_mass_sub = 340; % g, mass in unknown liquid
    iron_rho = 7.8; % g/cm^3
    fprintf("submerged iron has 3 forces...bouy, g, and balancer\n");
    fprintf("F_bouy = weight of fluid displaced by the iron && mass_air - mass_sub\n");
    fluid_mass = iron_mass_air - iron_mass_sub; % mass is conserved
    fprintf("a) fluid mass displaced =  %0.1f g   no decimal\n", fluid_mass);

    % find volume from density
    iron_volume = iron_mass_air / iron_rho;
    fprintf("b) V_iron =  %0.3f  g/cm^3\n", iron_volume);

    % for the unknown fluid
    fluid_rho = fluid_mass / iron_volume;
    fprintf("c) fluid_rho :  %0.3f  g/cm^3\n", fluid_rho);
end


%%%%%~~~~


if select == 9
    % laminar flow, pressure drop in a pipe
    p_drop = 95; % Pa
    radius = 7.5/1000; % m
    speed = 12/1000; %m/s
    area = pi * radius^2; % m^2

    % find net force
    fprintf("using: Delta(P) * area = F_net\n");
    f_net = p_drop * area;
    fprintf("a) F_net :  %0.5f  N\n", f_net);

    % find the power in mW
    fprintf("using power= force * velocity...\n");
    power = f_net * speed;
    fprintf("b) power:  %0.4f  mW\n", power*1000);

end


%%%%%~~~~


if select == 10
    Q = (175/60)/1000; % L/min --> L/s --> m^3/s 
    l = 45; % m
    d = 7.5/100; % m
    p_pump = 7.5e6; % N/m^2  "Pa"

    % from  Q = delat(P) / R   ...assume atmosphere same
    R = p_pump/Q;
    fprintf("a)  resistance in hose:  %0.3f  Ns/m^5\n", R);

    % viscosity assuming laminar flow...  looks wrong on website
    temp_a = p_pump / Q;
    temp_b = pi * (d/2)^4;
    temp_c = 8 * l;
    viscosity = (temp_a * temp_b) / temp_c;
    fprintf("b)  eta visc:  %0.3f  Ns/m^2    ???\n", viscosity);

    % power neglects height...it is same, pow=delat(P) * Q
    watts = p_pump * Q;
    fprintf("c)  power:  %0.1f W  no decmial\n", watts);
end


%%%%%~~~~


if select == 11
    fprintf("flow rate should be same no matter what cross-section looks like..conserve mass\n");
    fprintf("a) volume flow rate must be same everywhere in the pipe\n");
    fprintf("volume/time = area * velocity...   if dV/dt is same, smaller area has more speed\n");
    fprintf("b) A_1 * v_1 = A_2 * v_2 = ...    so, find where smallest area is\n");
    fprintf("c) solve  A_1 * v_1 = A_3 * v_3   -->    v_3/v_1 = A_1/A_3\n");
    fprintf("use bernolli, no potential...  P + (1/2) [rho] v^2 = const\n");
    fprintf("so section with lowest kentic energy has the highest pressure\n");
    fprintf("d) pick section with largest cross-section area\n");
end


%%%%%~~~~


if select == 12
    % you want to find flow rate in cm^3/s
    % gas engine, travels at 100 km/h
    % avereges 10 km/L
    fprintf("solve  Q=V/t  ...keep units right\n");
    % Q = speed / gas milage
    speed_kmph = 100; % km/h
    milage_kmpl = 10; % km/L
    Q_Lph = speed_kmph / milage_kmpl; % L/H
    cm3_per_liter = 1000; % cm^3 / L
    hours_per_second = 1/3600; % H / s
    Q_cmps = Q_Lph * cm3_per_liter * hours_per_second;
    fprintf("Q=  %0.2f  cm^3/s\n", Q_cmps);

    % kind of a seperate problem
    v = 2; % m/s   water moves in hose
    d = 1.6/100; % m   internal diameter of hose
    A = pi*(d/2)^2; % m^2
    Q_m3_p_s = A*v;
    liters_in_m3 = 1000;
    Q_L_s = Q_m3_p_s * liters_in_m3;
    fprintf("Q=  %0.3f  L/s\n", Q_L_s);

    v_noz = 15; % m/s   fluid velocity in nozzel
    d_noz = sqrt(Q_m3_p_s / (v_noz * pi)) * 2;
    fprintf("nozzel diamter is:  %0.3f cm\n", d_noz*100)
end


%%%%%~~~~


if select == 13
    h = 221; % m    Hoover Dam is this tall
    P_out = 1300e6; % MW  power output from dam
    depth = 150; % water is taken in at 150 meters
    Q_avg = 650; % m^3/s    average flow rate
    rho_water = 1000; % kg/m^3

    % power in the flow...   P = [rho] * g * h "depth taken in" * Q
    P_flow = rho_water * 9.81 * depth * Q_avg;
    fprintf("a) P_flow = %0.2f  E8 W\n", P_flow/1e8);

    % ratio of input???
    fprintf("b) try:  %0.3f   ???\n", P_flow/680e6);
end


%%%%%~~~~


if select == 14
    % airplane wing should get 1,000N of lift / m^2
    v_bot = 60; % m/s at take off    ground speed relative to bottom of wing
    rho_air = 1.29; % kg/m^3
    p_delta = 1000; % ideal lift...in Pa, or N/m^2

    % how fast must air-speed on top of wing be to get this ideal lift?
    % the height is considered the same... P_1 + (1/2)[rho_1][g][v_1]^2 + [rho_1][g][h]...h drops
    % just using kentic energy
    % P_1 - P_2 = (1/2) [rho_air] ([v_2]^2 - [v_1]^2)    solve for v_2
    % solve for v_2...
    temp = (2 * p_delta) / rho_air;
    v_top = sqrt(temp + v_bot^2);
    fprintf("a)  air speed on top of wing is:  %0.1f  m/s\n", v_top);

    v_cruise = 245; % m/s
    rho_hi = rho_air/4; % reduces 1/4
    temp = (2 * p_delta) / rho_hi;
    v_top = sqrt(temp + v_cruise^2);
    fprintf("b)  air speed on top of wing is:  %0.1f m/s   no decimal\n", v_top);
end


%%%%%~~~~


if select == 15
  % fill up pipe, at angle
  theta_deg = 25; % from horizontal...
  theta_rad = deg2rad(theta_deg); 
  l = 0.55; % m
  Atank_Anoz = 150;  % nozzel cross section area is 150x smaller than tank...
  rho_water = 1000;

  % find the pressure, max absolute, at tip of gun, given it is angled and being loaded
  % pressure at depth, delta(P) = [rho][g][depth below surface]
  depth = l * sin(theta_rad);    % how deep the gun goes below surface
  p_deep = rho_water * 9.81 * depth;
  fprintf("a)  p_absolute = p_0 + p_deep =  p_0 + [rho_w][g][l]sin(theta_rad)\n");
  p_0 = 1.013e5; % Pa at sea level...
  p_abs = p_0 + p_deep;
  fprintf("b)  p_abs =  %0.3f Pa\n", p_abs);

  % pluger pulled, and filled
  p_tank = 0.11 * p_0;  % Pa
  fill_time = 41;       % s
  % P1 + (1/2)[rho][v1]^2 + [rho][g][y1] = P2 + (1/2)[rho][v2]^2 + [rho][g][y2]
  % 1== nozzel, 2== tank
  % P1 = P2 + (1/2)[rho]{[v2]^2 - [v1]^2} + [rho][g][y2 - y1]
  v_tank = l / fill_time;  % fills gun evenly...  length / time 
  % y2 - y1 = l sin(theta)   ...  depth
  v_noz = v_tank * Atank_Anoz; % use the given ratio, solve   A_tank * v_tank = A_noz * v_noz
  p_noz = p_tank + (1/2) * rho_water * (v_tank^2 - v_noz^2) + rho_water * 9.81 * depth;
  fprintf("c) pressure at nozzel:  %0.3f Pa\n", p_noz);
  
end


%%%%%~~~~


if select == 99
    fprintf("\n\tDONE\n");
end


%%%%%%%%~~~~~~~~END>  v1_ch14_flud.m

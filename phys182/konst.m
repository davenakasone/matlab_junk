%{
    fyi...

    V_sphere= (4/3)*pi*r^3
    A_circle = pi*r^2

%}

classdef konst
    properties (Constant)
        c = 299.8e6;                       % speed of light                                     [m/s]
        boltzmann = 1.38e-23;              % J/K  Boltzmann's constant                          [?]
        g = 9.8;                           % gravitational accelaration at surface of earth     [m/s^2]
        cm_per_inch = 2.54;                % there are 2.54 cm in one inch                      [cm/in]
        pound_2_newton = 4.448;            % one pound is 4.448 Newtons                         [N/lb]
        newton_2_pound = 1/4.448;          % one Newton is (1/4.448) pounds of force            [lb/N]
        rho_water = 1000;                  % density of water                                   [kg/m^3]
        pa_2_atm = 1/1.013e5;              % 1 Pa to atmospheres                                [atm/Pa]
        atm_2_pa = 1.013e5;                % 1 atmosphere to Pa                                 [Pa/atm]
    end
end
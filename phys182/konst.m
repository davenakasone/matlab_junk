%{
    fyi...
%}

classdef konst
    properties (Constant)
        avogadro = 6.02e23;                % molecules in a mole                                [1/molecules]
        c = 299.8e6;                       % speed of light                                     [m/s]
        boltzmann = 1.38e-23;              % J/K  Boltzmann's constant                          [J/K]
        boltzmann_radiation = 5.67e-8;     % Boltzmann for radiation                            [J/(s m^2)]
        g = 9.8;                           % gravitational accelaration at surface of earth     [m/s^2]
        cm_per_inch = 2.54;                % there are 2.54 cm in one inch                      [cm/in]
        pound_2_newton = 4.448;            % one pound is 4.448 Newtons                         [N/lb]
        newton_2_pound = 1/4.448;          % one Newton is (1/4.448) pounds of force            [lb/N]
        rho_water = 1000;                  % density of water                                   [kg/m^3]
        pa_2_atm = 1/1.01e5;              % 1 Pa to atmospheres                                [atm/Pa]
        atm_2_pa = 1.01e5;                % 1 atmosphere to Pa                                 [Pa/atm]
        kcal_2_joule = 4186;               % 1 kcal to 4186 J                                   [J/kcal]
        joule_2_kcal = 1 / 4186;           % 1 kcal is (1/4186) J                               [kcal/J] 
        btu_2_joule = 1055.1;              % 1 BTU has 1055.1 J, 1lb water 1 F\deg              [J/btu]
        joule_2_btu = 1/1055.1;            % 1 J is 1/1055.1 BTU                                [but/J]
        gas_const_J_mol_K = 8.31;          % ideal gas constant                                 [J/(mol K)]
        gas_const_cal_mol_K = 1.99;        % ideal gas constant                                 [cal/(mol K)]
        gas_const_L_atm_mol_K = 0.821;     % ideal gas constant                                 [(L atm)/(mol K)]
    end
    methods (Static)
        function a_circle = area_of_circle(radius)
            a_circle = pi * radius^2;
        end
        function v_sphere = volume_of_sphere(radius)
            v_sphere = (4/3) * pi * radius^3;
        end
        function fahrenheit_out = celsius_2_fahrenheit(celsius_in)
            fahrenheit_out = (9/5) * celsius_in + 32;
        end
        function celsius_out = fahrenheit_2_celsius(fahrenheit_in)
            celsius_out = (5/9) * (fahrenheit_in - 32);
        end
        function kelvin_out = celsius_2_kelvin(celsius_in)
            kelvin_out = celsius_in + 273; % 0.15...
        end
        function celsius_out = kelvin_2_celsius(kelvin_in)
            celsius_out = kelvin_in - 273; % 0.15...
        end
        function kelvin_out = fahrenheit_2_kelvin(fahrenheit_in)
            kelvin_out = (5/9) * (fahrenheit_in - 32) + 273.15;
        end
        function fahrenheit_out = kelvin_2_farhrenheit(kelvin_in)
            fahrenheit_out = (9/5) * (kelvin_in - 273) + 32; % 0.15...
        end
        function speed = max_bolt_v(v, m, T)
            temp_a = 4 / sqrt(pi);
            temp_b = 2 * (1.38e-23) * T;
            temp_c = (m / temp_b)^(3/2);
            temp_d = (v^2) * exp( (-1 * m * v^2) / temp_b);
            speed = temp_a * temp_c * temp_d;
        end
    end
end
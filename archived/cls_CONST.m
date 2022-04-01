classdef cls_CONST
    % keep common constants here
    % make an obj instance, reference consts at will + basic conversions
    
    properties
    end
    
    properties (Constant)
        
        % prefixes
        yotta = 10^24;  % Y yotta
        zetta = 10^21;  % Z zetta
        exa   = 10^18;  % E exa
        peta  = 10^15;  % P peta
        tera  = 10^12;  % T tera
        giga  = 10^9;   % G giga
        mega  = 10^6;   % M mega
        kilo  = 10^3;   % k kilo
        hecto = 10^2;   % h hecto
        deka  = 10^1;   % da deka 
        deci  = 10^-1;  % d deci
        centi = 10^-2;  % c centi
        milli = 10^-3;  % m milli
        micro = 10^-6;  % u micro
        nano  = 10^-9;  % n nano
        pico  = 10^-12; % p pico
        femto = 10^-15; % f femto
        atto  = 10^-18; % a atto 
        zepto = 10^-21; % z zepto
        yokto = 10^-24; % y yokto
        
        avo = 6.022e23;                  % Avogadro's number     atoms per  mol, any element
        
        % break-down strength E in V/m
        % only valid at low freq < 1 kHz   may change otherwise
        bd_bat = 7.5e6;         % Barium Titanate
        bd_pap = 12e6;          % Paper
        bd_gls = 35e6;          % Glass 
        bd_mca = 70e6;          % Mica
        bd_bkl = 20e6;          % Bakelite
        bd_qtz = 30e6;          % Quartz (fused)
        bd_rub = 25e6;          % Rubber(hard)
        bd_paf = 30e6;          % Paraffin
        bd_pet = 12e6;          % Petrolium Oil
        bd_air = 3e6;           % Air  @ 1 ATM
        
        bolt = 8.62e-5;   % boltzmann constant in electron volts (1.6e-19) / degree kelvin
        
        c = 3e8;                         % speed of light  m/s
        
        %conductivity is siemens/meter  <->  1 / ohm m    @  20°C
        cond_slv = 6.1e7;                % Silver
        cond_cop = 5.8e7;                % Copper
        cond_gld = 4.1e7;                % Gold
        cond_alm = 3.5e7;                % Aluminum
        cond_tng = 1.8e7;                % Tungsten
        cond_znc = 1.7e7;                % Zinc
        cond_brs = 1.1e7;                % Brass
        cond_irn = 10^7;                 % Iron 
        cond_led = 5e6;                  % Lead
        cond_mrc = 10^6;                 % Mercury
        cond_car = 3e4;                  % Carbon
        cond_swa = 4;                    % Sea Water
        cond_grm = 2.2;                  % Germanium
        cond_sil = 4.4e-4;               % Silicon
        cond_dwa = 10^-4;                % Distilled Water
        cond_eth = 10^-5;                % Dry Earth
        cond_bkl = 10^-10;               % Bakelite
        cond_pap = 10^-11;               % Paper
        cond_gls = 10^-12;               % Glass
        cond_por = 10^-12;               % Porcelin
        cond_mca = 10^-15;               % Mica
        cond_paf = 10^-15;               % Paraffin 
        cond_rub = 10^-15;               % Rubber (hard)
        cond_qtz = 10^-17;               % Quartz (fused)
        cond_wax = 10^-17;               % Wax
        
        % ep_r: relative permittivity aka "dielectric const"  no units
        % only valid at low freq < 1 kHz   may change otherwise
        epr_bat = 1200;         % Barium Titanate
        epr_swa = 80;           % Sea Water
        epr_dwa = 81;           % Distilled Water
        epr_nyl = 8;            % Nylon
        epr_pap = 7;            % Paper
        epr_gls_min = 5;        % Glass (minimum)
        epr_gls_max = 10;       % Glass (maximum)
        epr_mca = 6;            % Mica
        epr_por = 6;            % Porcelain
        epr_bkl = 5;            % Bakelite
        epr_qtz = 5;            % Quartz (fused)
        epr_rub = 3.1;          % Rubber(hard)
        epr_sox = 3.9;          % Silicon Dioxide
        epr_wod_min = 2.5;      % Wood(minimum)
        epr_wod_max = 8;        % Wood(maximum)
        epr_pst = 2.55;         % Polystyrene 
        epr_ppp = 2.25;         % Polypropylene
        epr_paf = 2.2;          % Paraffin
        epr_pet = 2.1;          % Petrolium Oil
        epr_air = 1;            % Air  @ 1 ATM
        epr_cop = 1;            % Copper
        
        ef = 1.6e-19;                    % fundamental unit of charge, coulumbs...  * -1 for electron
        
        ep0 = (10^(-9))/(36*sym(pi));    % permittivity of free space (really 8.854e-12)  F/m

        g = 9.8;                         % accelaration due to force of gravity ,earth's surface m/s^2
        
        G = 6.67e-11;                    % universal gravitational constant  N (m/kg)^2

        k = 9e9;                         % proportionality constant  1 / ( 4 pi ep0 )   m/F
        
        mu0 = 4*sym(pi)*(10^-7);         % permeability of free space H/m
        
        % silicon 2, 8, 4(valance), 14 total electrons
        ep_sox = 3.45e-11;               % permitivitty of silicon dioxide F/m
        sil_adense = 5e22;               % atomic density of silicon   atoms/cm^3
        sil_amu = 28;                    % atomic weight...28 nuetrons / 28 protons
        silB = 7.3e15;                   % paramter of silicon  use for temperature cm^-3 K^-3/2
        sil_dense = 2.33;                % density of silicone g/cm^3
        sil_Dp = 12;                     % hole difusion const   cm^2 / s
        sil_Dn = 35;                     % electron difusion const  cm^2 / s
        sil_Eg = 1.12;                   % band gap energy of silicon  eV   ...1.6e-19
        sil_mu_p = 480;                  % hole molbility const, intrinsic silicon  cm^2 / V s
        sil_mu_n = 1350;                 % electron mobility const, instrinic silicon  cm^2 / V s
        sil_ni_300 = 1.5e10;             % silicon  intrinsic carrier density  T = 300 K , carriers / cm^3
        sil_perm = 1.04e-12;             % permitivity of silicon...2 * ep0    F / cm
        sil_rho = 5e22;                  % silicon atoms per unit volume   atoms/cm^3
        sil_VTroom = .025;               % short-hand V_T for silicon at 20°C     V @ 20°C...25 mV
        VT300 = 25.9e-3;                 % thermal voltage for Eisten relationship room temp 20C°
        
        % ur : relative permeability 
            % diamagnetic
            ur_bis = .999833;    % Bismuth
            ur_mer = .999968;    % Mercury
            ur_sil = .9999736;   % Silver
            ur_led = .9999831;   % Lead
            ur_cop = .9999906;   % Copper
            ur_wat = .9999912;   % Water
            ur_hyd = 1.0;        % Hydrogen(STP)
            % paramagnetic
            ur_oxy = .999998;    % Oxygen(STP)
            ur_air = 1.00000037; % Air
            ur_alm = 1.000021;   % Aluminum
            ur_tng = 1.00008;    % Tungsten
            ur_plt = 1.0003;     % Platinum
            ur_mng = 1.001;      % Manganese
            % ferromagnetic
            ur_cbt = 250;        % Cobalt
            ur_nkl = 600;        % Nickel
            ur_irn = 5000;       % Iron (soft)
            ur_sfe = 7000;       % Silicon-Iron 
    end
    
    properties (Access = protected)
        obj_name;
    end
    
    methods
        
        function this_obj = cls_CONST()    %  constructor
            this_obj.obj_name = "const";
            %fprintf('\t use the obj for consts\n');
            %fprintf('see const property list and make sure you have updated consts\n');
            %fprintf('be careful for unit changes / anything not SI\n');
        end
        
        function check(this_obj)
            fprintf('\n\t%s   can still be found in memory\n', this_obj.obj_name);
        end
        
        function help(this)
            fprintf('this CONST object has numeric values and corresponsing unit values as strings\n');
            fprintf('this object currently has:\n');
            fprintf('the speed of light is { consts.ls } = %d   { consts.ls_u } %s\n', this.ls, this.ls_u);
        end
        
        function ni = ni4sil(this, kelv)
            ni = this.silB * kelv^(3/2) * exp( -this.sil_Eg / ( 2 * this.bolt * kelv ) );
        end
  
    end 
    
    
    methods(Static)
        
        function kelv = far2kel(farIn)
            cel = (5/9)*farIn - (160/9);
            kelv = cel + 273; % maybe change to 273.15
        end
        
        function kelv = cel2kel(celIn)
            kelv = celIn + 273; % maybe change to 273.15
        end
        
        function vt = kel2vt(kelIn)
            vt = (.0862 * kelIn)/1000; % short hand from p.215, ch4   in volts   
            % use 25 mV if room temeperature ( 20C ) " sil_vtRoom"
        end
        
        function vt = cel2vt(celIn)
            vt = ( .0862 * (273+celIn) )/1000; % short hand from p.215, ch4   in volts   
            % use 25 mV if room temeperature ( 20C ) " sil_vtRoom"
        end
        
        function par_R = para(R1, R2)
            par_R = (R1 * R2) / (R1 + R2);
        end
        
        function alp_out = alp_cal(om, mu, ep, sg) % attenuation coef of EM wave[planar TEM]
            t1 = (mu*ep)/2;
            t2 = (sg/(om*ep))^2;
            alp_out = om * sqrt( t1 * (sqrt(1+t2) - 1)  ); % lossy...not lossless
        end
        
        function bet_out = bet_cal(om, mu, ep, sg) % wave number EM wave[planar TEM]
            t1 = (mu*ep)/2;
            t2 = (sg/(om*ep))^2;
            bet_out = om * sqrt( t1 * (sqrt(1+t2) + 1)  );  % lossy, general case
        end
        
        function eta_out = eta_cal(om, mu, ep, sg) % intrinsic impedance of EM wave
            eta_out = sqrt( (1j*om*mu)/(sg + 1j*om*ep) );
        end
        
    end
    
end


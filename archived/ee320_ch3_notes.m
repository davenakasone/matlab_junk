%{
    ee320_ch3_notes
    
    #1      ex3.1    n_i  , intrinic carrier density of silicon   room temp 300k
    #2      pp3.1    n_i  , intrinic carrier density of silicon    50k, 350k   big change
    #3      ex3.2 , pp3.2   n-type silicon   Nd --> nn, pn
    #4      pp3.3   boron doped silicon
    #5      ex3.3   resistivity  regular silcon and ptype
    #6      pp3.4   resistivity, current, ect
    #7      ex3.4   difusion current
    #8      pp3.5   difusion to obj  use graph
    #9      ex3.6   eistein const
    #10     ex3.5   pn junciton at room temperature
    #11     ex3.6   pn junciton on voltage source

    #991    hw3.1
    #992    hw3.2
    #996    hw3.6
    #998    hw3.8



%}
%color = uisetcolor([1, 1, 0], 'Selecf Color');    % [.9, .9, .9] is a nice gray
%clear all;  if things aren't going away
clc;
close all;
clearvars;

                sel = 998;  % CHANGE CHANGE CHANGE


% EE330 globals
global rx; global ry; global rz; % rectangular params  rx, ry, rz        
global cr; global cf; global cz; % cylindrical params  cr, cf, cz
global sr; global st; global sf; % spherical params    sr, st, sf
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); syms cf; assume(cf,'real'); syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); syms st; assume(st,'real'); syms sf; assume(sf, 'real');
global pL; syms pL; assume(pL, 'real'); % line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); % surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); % volume charge denisty in C / m^3

% Math432 globals   only X и Y are needed in some helper funcitons 
global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Zxy; syms Zxy; Zxy = X + 1j*Y;% just use it as a hold and imply  when needed
global Z; syms Z;            % your temporary Z for integration and differentiation
global Ur; syms Ur; assume(Ur, 'real');   % Ur ( X, Y)  ...real part of Z = Ur(X,Y) + j Vi(X,Y)      
global Vi; syms Vi; assume(Vi, 'real');   % Vi ( X, Y)  ...imag part of Z = Ur(X,Y) + j Vi(X,Y)     
% your function = Ur + j Vi   OR your function = operation(Z)
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var
global pu; syms pu; assume(pu, 'real'); % paramater u for surface trace
global pv; syms pv; assume(pv, 'real'); % paramater v for surface trace
% series 
global in; syms in; assume(in, {'real', 'integer'}); % index n
%global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global ik; syms ik; assume(ik, {'real', 'integer'}); % index k, usually for convolution
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
global T; syms T; assume(T, 'real');  % period bound, continous case
global t0; syms t0; assume(t0, 'real'); % arbitrary t0, usually offset
global tau; syms tau; assume(tau, 'real'); % dummy tau, usually for convolution

ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #1
if sel == 1
    temp = 300; % degrees K
    ni = const.Bsil * temp^(3/2) * exp( -const.Eg_sil / ( 2 * const.bolt * temp ) );
    fprintf(' silicon is about %.1f * 10^10  carriers / cm^3\n', ni *(10^-10));
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    K = 50;
    ni = const.Bsil * K^(3/2) * exp( -const.Eg_sil / ( 2 * const.bolt * K ) );
    fprintf('intrinsic carrier density ni, T = %.0f K : %d\n', K, ni );
    K = 350;
    ni = const.Bsil * K^(3/2) * exp( -const.Eg_sil / ( 2 * const.bolt * K ) );
    fprintf('intrinsic carrier density ni, T = %.0f K : %d\n', K, ni );
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    Nd = 10^17;  % donor atoms / cm^3
    k = 300; % temperature is 300 kelvin
    nn = Nd;  % free electrons are about = impurity concentration
    pn = ( const.ni_sil_300 )^2 / Nd;  % holes are ni^2 / doping concentration
    fprintf('n-type holes = %d    holes / cm^3  notice nn >> ni,  nn >> pn\n', pn);
    ni = const.ni4sil(350);
    pn = (ni^2)/Nd;
    fprintf('n-type holes = %d    holes / cm^3  notice nn >> ni,  nn >> pn\n', pn);
    fprintf(' no change in n_n , but a lot more holes\n');
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    k = 300;
    ni = const.ni_sil_300;
    np = ni/10^6; % n_p drops by a large factor
    pp = (ni^2)/np;
    fprintf('with drop in np, pp must be %s   atoms/cm^3\n', pp);
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    p = const.ni_sil_300; % p = n = ni for intrinsic silicon
    rho  = 1 / ( const.ef * ( p * const.mu_p_sil + p * const.mu_n_sil) );
    fprintf('the resistivity of intrinic silicon is %d  ohm * cm\n', rho);
    Na = 10^16; % acceptors / cm^3
    pp = Na;
    np = ( const.ni_sil_300 )^2 / Na;
    rho = 1 / ( const.ef * ( pp * 400 + np * 1110));  % mu_p  mu_n given
    fprintf('p-type silicon resistivity = %.3f  ohm * cm\n', rho);
    fprintf('doping reduced resistivity a lot and it is mostly determined by concentration\n');
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % n type silicon bar, L = 2e-6 m   1V across  
    Nd = 10^16;
    mu_n = 1350;
    V = 1;
    L = (2e-6)*100; % length in centimeters
    E = V/L; % volts per cm
    nn = Nd; % for n-type
    pn = (const.ni_sil_300^2)/Nd;
    vn_drift = -mu_n*E;
    fprintf('electron drift velocity (against E) = %d   cm / s\n', -1*vn_drift);
    time = L / vn_drift;
    fprintf('electron crosses length in %.2f  ps\n', -1*time*(1e12));
    Jsn = -1*const.ef * vn_drift * nn;
    fprintf('drift current density Jsn = %d   A / cm^2\n', Jsn);
    Am = (.25e-6); % given cross section area = .25um^2, convert to cm^2
    Acm = Am*(10^4);%(sqrt(Am)*(1e2))^2;
    I = Jsn * Acm
    rho = E / Jsn;
    R = rho*L/Acm % something is off
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    p0 = 10^16; % holes / cm^3
    Lp = (1e-6)*(10^2); % length to cm
    Dp = 12; % cm^2 / s
    A = (100e-6)*(1e-2); % um^2 --> cm^2
    px = p0*exp(-rx/Lp);
    dpx_dx = diff(px, rx, 1);
    dif = subs(dpx_dx, rx, 0);
    Jdp = -const.ef * Dp * dif; 
    fprintf('diffusion current density Jdp %d  A/cm^2\n', Jdp);
    Idp = Jdp * A;
    fprintf('diffusion current Idp %.3f  uA\n', Idp*(1e6));
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    n0 = (10^17)*(1e-6); % free electrons / cm^3  --> m^3
    W = (.5e-6); % um --> m^3
    dnx_dx = ((0-n0) / (W-0)); %*(1e12); % --> m^2   
    Jdn = ( const.ef * (1e-4)*const.sil_Dn * dnx_dx ) * (1e6) % A/m^2    should be 112uA/m^2
    req_I = .001; % 1mA needed
    A = req_I/(112e-6);
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    Dn = const.sil_mu_n * const.VT300;    % by the eistein relationp
    Dp = const.sil_mu_p * const.VT300;
    fprintf('Dp must be %.3f  and Dn must be %.3f    cm^2/s\n', Dp, Dn);
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    NA = 10^18; % holes / cm^3
    ND = 10^16; % free electrons / cm^3
    A = 10^-4;  % cm^2   cross section area
    ni = 1.5e10; % given
    pp = NA;
    fprintf('holes pp should be about equal to acceptor concentration NA: %d\n', pp);
    nn = ND;
    fprintf('free electrons nn should be about donor concentration    ND: %d\n', nn);
    np0 = (ni^2)/NA; % or ni^2 / pp
    fprintf('minority free electrons on p-side: %.3f / cm^3\n', np0); 
    pn0 = (ni^2)/ND; % or ni^2 / nn
    fprintf('minority holes on n-side: %.3f / cm^3\n', pn0);
    V0 = const.VT300 * log(NA*ND/ni^2);
    fprintf('voltage across depletion zone = %.3f V\n', V0);
    W = sqrt((2*V0*const.sil_perm/const.ef)*((1/NA)+(1/ND)));
    fprintf('width is %.6f cm   or %.3f  um\n', W, W*(1e6));
    xn = W*(NA/(NA+ND));
    xp = W*(ND/(NA+ND));
    fprintf('width goes xp = %.3f  um   to xn = %.3f um  for pn junction\n', xp*(1e4), xn*(1e4));
    QJ = A*xn*ND*const.ef; % lots of ways to get it
    QJ = A*xp*NA*const.ef;
    fprintf('charge on either side: %.3f pC\n', QJ*(1e12));
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    NA = 10^18; % holes / cm^3 
    ND = 10^16; % free electrons / cm^3 
    A = 10^-4; % cm^2  --> m^2
    ni = 1.5e10; % carriers / cm^3  
    Lp = (5e-6)*100; % m --> cm             be careful with these guys -->cm for all
    Ln = (10e-6)*100; % m --> cm
    Dp = 10; %cm^2 / Vs      n-region
    Dn = 18; % "        "    p-region
    I = .1e-3; % A
    Is = A * const.ef * (ni^2) * ( (Dp/(Lp*ND)) + (Dn/(Ln*NA)) );
    % if I = Is(exp(V/Vt) - 1).... about Is exp(V/VT)
    V = const.VT300 * log(I/Is);
    Ip = A * const.ef * (Dp/Lp)*(ni^2/ND)*(exp(V/const.VT300)-1);
    In = A * const.ef * (Dn/Ln)*(ni^2/NA)*(exp(V/const.VT300)-1);
    ratio = Ip /In; % 111x stronger  ... doping was 111x in p-type compared to n
end


%------------------------------------------------------------------------------------------ #991
if sel == 991
    B       = 7.3e15;                                %   1 / cm^3 Kelvin^(-3/2)
    Eg      = 1.12;                                  %   electron volts
    k       = 8.62e-5;                               %   electron volts / Kelvin  "Boltzman const"
    dense   = 5e22;                                  % silicon atoms / cm^3
    celcius = [-55, 0, 20, 75, 125];                 % degrees celcius
    ni      = [0  , 0, 0 , 0 , 0];
    kelvin  = [0  , 0, 0 , 0 , 0];
    frac    = [0  , 0, 0 , 0 , 0];
   
    for idx = 1:5
        kelvin(idx)  = celcius(idx) + 273; % degrees Kelvin
        % carriers / cm^3
        ni(idx)      = B * kelvin(idx)^(3/2) * exp( -Eg/ ( 2 * k * kelvin(idx)) );
        frac(idx)    = ( ni(idx) / dense );            % fraction of atoms ionized 
        fprintf(' at %.1f °C { %.1f K } :\n', celcius(idx), kelvin(idx));
        fprintf(' \t ni = %d  carriers / cm^3\n', ni(idx))
        fprintf(' \t\t %d  fraction of atoms ionized to non-ionized\n\n', frac(idx));
    end
    result = zeros(4, 5);
    result(1,:) = celcius;
    result(2,:) = kelvin;
    result(3,:) = ni;
    result(4,:) = frac;
    figure('Position', [20, 20, 900, 900]);
    hold on;
    xlabel('K', 'FontSize', 14);
    ylabel('ni [ carriers / cm^3 ]', 'FontSize', 14);
    title(' temperature (K) vs carriers / cm^3', 'FontSize', 16);
    plot(kelvin, ni, 'r--', 'LineWidth', 2);
    for idx = 1:5
        plot(kelvin(idx), ni(idx), 'b.', 'MarkerSize', 20);
    end
    hold off;
    figure('Position', [20, 20, 900, 900]);
    hold on;
    xlabel('K', 'FontSize', 14);
    ylabel('fraction of ionized atoms', 'FontSize', 14);
    title(' temperature (K) vs fraction of ionized atoms', 'FontSize', 16);
    plot(kelvin, frac, 'g--', 'LineWidth', 2);
    for idx = 1:5
        plot(kelvin(idx), frac(idx), 'm.', 'MarkerSize', 20);
    end
    hold off;
end


%------------------------------------------------------------------------------------------ #992
if sel == 992
    ni_Si = 1.5e10;      % given
    kelvin = 300;        % temperature
    B      = 3.56e14;    % 1 / cm^3 K^(3/2)
    Eg     = 1.42;       % electron volts
    k      = 8.62e-5;    % electron volts / Kelvin
    ni_GaAs = B * kelvin^(3/2) * exp(-Eg / (2 * k * kelvin));
    fprintf('silicon has %s  carriers / cm^3\n', ni_Si);
    fprintf('gallium arsenide has %d carriers / cm^3\n', ni_GaAs);
    fprintf(' \tsilicon has %d times the carriers of gallium arsenide\n', ni_Si/ni_GaAs);
end


%------------------------------------------------------------------------------------------ #996
if sel == 996
    B   = 7.3e15;      %   1 / cm^3 Kelvin^(-3/2)
    Eg  = 1.12;        %   electron volts
    k   = 8.62e-5;     %   electron volts / Kelvin  "Boltzman const"
    NA  = 10^17;       % boron doped silicon , holes / cm^3
    
    kelvin = 27+273; 
    ni = B * kelvin^(3/2) * exp(-Eg / ( 2 * k * kelvin ) ); % about 1.5e10 as given
    pp = NA;
    np = (ni^2)/NA;
    fprintf('at 300 K, hole concentration pp = %s holes / cm^3\n', pp);
    fprintf('\t electron concentratin np = %d  free electrons / cm^3\n\n', np);
    
    kelvin = 125+273; 
    ni = B * kelvin^(3/2) * exp(-Eg / ( 2 * k * kelvin ) ); % about 1.5e10 as given
    pp = NA;
    np = (ni^2)/NA;
    fprintf('at 398 K, hole concentration pp = %s holes / cm^3\n', pp);
    fprintf('\t electron concentratin np = %d  free electrons / cm^3\n\n', np);
end


%------------------------------------------------------------------------------------------ #998
if sel == 998
    un = 1350 * (10^-4); % cm^2 / V s  --> m^2 / V s
    up = 480 * (10^-4);  % cm^2 / V s     --> m^2 / V s
    L = 10e-6; % length in meters
    V = 3; % volts
    E = (V/L); % electric field, V/m
    vp_drift = up * E; % m/s
    vn_drift = -un * E; % m/s
    fprintf('hole drift velocity is     %d  m / s\n', vp_drift);
    fprintf('electron drift velocity is %d m / s\n', vn_drift);
    fprintf('\t\t the electron moves faster, but in oppisite direction\n');
end
    
    
    
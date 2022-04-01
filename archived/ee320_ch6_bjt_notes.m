%{
    ee320_ch6_bjt_notes

    #1      pp6.1   infer Vbe
    #2      pp6.2   work with alp/bet
    #3      pp6.3   the 3 currents
    #4      pp6.4   npn with some parameters
    #5      pp6.5   more params
    #6      pp6.6   working ex1
    #7      pp6.7   compare saturation currents
    #10     pp6.10  pnp
    #11     p6.11   pnp

    
%}
clc;
close all;
clearvars;
sympref('MatrixWithSquareBrackets', 1);
sympref('PolynomialDisplayStyle', 'ascend');


                sel = 11;  % CHANGE CHANGE CHANGE

                
% EE320 globals
global alpha; syms alpha; assume(alpha, 'real');
global beta; syms beta; assume(beta, 'real');

% EE330 globals
global sig; syms sig; assume(sig, 'real'); % resistivity , sub in and out with const
global c1; syms c1; % just a holder for integration constant, solve for latter
global c2; syms c2; % just a holder for integration constant, solve for latter
global c3; syms c3; % just a holder for integration constant, solve for latter
global c4; syms c4; % just a holder for integration constant, solve for latter
global rx; global ry; global rz; % "rectangular x"   "rectangular y"  "rectangular z"       
global cr; global cf; global cz; % "cylindrical rho" "cylindrical fi" "cylindrical z"       fi = phi 
global sr; global st; global sf; % "spherical radius" "spherical theta" "spherical fi"   to avoid confusion
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); assume(cr >= 0);
syms cf; assume(cf, 'real'); assume(cf >= 0);
syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); assume(sr >= 0);
syms st; assume(st, 'real'); assume(st >= 0);
syms sf; assume(sf, 'real'); assume(sf >= 0);
global pL; syms pL; assume(pL, 'real'); assume(pL >= 0); % line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); assume(pS >= 0); % surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); assume(pV >= 0); % volume charge denisty in C / m^3
global ep; syms ep; assume(ep, 'real'); % as in ep = ep0 * epr    or D = ep * E
global ep1; syms ep1; assume(ep1, 'real'); % represent permittivity of first region
global ep2; syms ep2; assume(ep2, 'real'); % represent permittivity of second region
global ep3; syms ep3; assume(ep3, 'real'); % if you need more than 3, then someone hates you
global epr; syms epr; assume(epr, 'real'); % as in epr = ep/ep0  or   epr = 1 + chie
global ep0; syms ep0; assume(ep0, 'real'); % sub out with const.ep0 at end of calculation (if needed)
global chie; syms chie; assume(chie, 'real'); % as in  chie = epr - 1   or chie = ep/ep0 - 1

% EE360 / Math432 globals   only X и Y are needed in some helper funcitons 
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
global s; syms s; % transform holder " j w0 "
global freq; syms freq; assume(freq, 'real'); % in Hz  f = 1 / T   or N
global omg; syms omg; assume(omg, 'real'); % omega , angular freq ... 2 pi f
global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global ik; syms ik; assume(ik, {'real', 'integer'}); % index k, usually for convolution
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
global T; syms T; assume(T, 'real');  % period bound, continous case
global t0; syms t0; assume(t0, 'real'); % arbitrary t0, usually offset
global tau; syms tau; assume(tau, 'real'); % dummy tau, usually for convolution

ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #0
if sel == 0
end


%------------------------------------------------------------------------------------------ #1
if sel == 1
    Vbe = .7;
    ic = 1e-3;
    ica = .1e-3;
    icb = 10e-3;
    
    Is = ic / ( exp( Vbe/const.sil_vtRoom) );
    Vbe1 = Vbe + const.sil_vtRoom * log( (ica/Is)/(ic/Is) ); % 0.6424 V
    Vbe2 = Vbe + const.sil_vtRoom * log( (icb/Is)/(ic/Is) ); % 0.7576 V
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    b1 = 50;
    b2 = 150;
    a1 = b1/(b1+1); % 0.9804
    a2 = b2/(b2+1); % 0.9934
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    ib = 14.46e-6;
    ie = 1.46e-3;
    Vbe = .7;
    
    ic = ie - ib; % 0.0014
    Is = ic / exp( Vbe / const.sil_vtRoom ); % 9.9950e-16
    b = ic/ib;  % 99.9682
    a = b/(b+1);  %  0.9901
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    a1 = .99;
    a2 = .98;
    b1 = a1/(1-a1); % 99.0000
    b2 = a2/(1-a2); % 49.0000
    ic = 10e-3;
    ib1 = ic/b1; % 1.0101e-04
    ib2 = ic/b2; % 2.0408e-04
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    Is = 10^-16;
    b = 100;
    a = b/(b+1); % 0.9901
    ic = 1e-3;
    Vbe = const.sil_vtRoom * log(ic/Is); % 0.7483 V
    Ise = Is/a; % 1.0100e-16  A
    Isb = Is/b; % 1.0000e-18  A
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % just remember Vce > Vbe to stay active
    R = (5-.69)/(1e-3);  % 4310 ohms 
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    Is = 10^-15;
    a = 1;
    Ae = 100;
    Isc = Is * Ae;  % 1.0000e-13
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    b = 50;
    a = b/(b+1);
    Is = 10^-14;
    ie = 2e-3;
    Veb = const.sil_vtRoom * log(a * ie / Is ); % 0.6500
    ic = Is * exp(Veb/const.sil_vtRoom) %   0.0020
    ib = ie - ic; % 3.9216e-05
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    Is = 10^-11;
    b = 100;
    a = b/(1+b);
    ic = 1.5;
    Veb = log(ic/Is)*const.sil_vtRoom; %  0.6433 V
end
    
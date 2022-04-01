%{
    ch6, electrostatic boundary-value problems    notes

    #1      pp6.1       Poisson
    #2      pp6.3       see the plates, use Laplace, put it cylindrical
    #3      pp6.4       this is a cone, keep it all on st, and Laplace it  

    #333    hw5, question 2

%}
clc;
close all;
clearvars;


                sel = 333;  % CHANGE CHANGE CHANGE


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

ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #1
if sel == 1
    syms rho_0; assume(rho_0, 'real');
    syms a; assume(a, 'real');
    
    pV = (rho_0 * rx) / a; % given volume charge density   ...it's a simple one-D device
    % bc1:  E = [0,0,0] @ x = 0
    % bc2: V = 0  @ x = a
    
    % all the action happens [1,0,0] normal
    % since rho_v is not 0 (there is free charge)  use   (d/dx)(dV/dx) = rho_v/ep
    int1 = int( -pV/ep, rx);
    int1 = int1 + c1;
    int2 = int(int1, rx);
    int2 = int2 + c2;
    E = ee.getGradRec(int2);
    c1f = 0; % solved
    c2f = (rho_0 * a^2) / (6 * ep); % solved
    V = subs(int2, [c1, c2], [c1f, c2f] );
    E = subs(E, [c1, c2], [c1f, c2f] );
        pretty(V);
        pretty(E); 
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    epr = 1.5; % between
    V = (200/sym(pi))*cf; % easy enough to see and fill in with given bondary conditions... all cf
    E = -1 * ee.getGradCyn(V); % [0, -200/(cr*pi), 0]
    D = (const.ep0 * epr) .* E; % [0, -300/(cr*pi), 0]
    % use the fact rho_s = D
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >  VECTOR
    % dS = [ 0, 1, 0]
    q = double(int(int( D(2), cr, .004, 1), cz, 0, 5))*(1e9) 
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % V(st = pi/4) = 50        V(st=pi/2) = 0     ...tiny gap, cone is charged, z=0 plane ground
    % this can be done completley spherical...which is good
    pt_rec = [-3,4,2];
    ee.feed(pt_rec, 'R');
    pt_sph = ee.pts(3,:);
    sta = pi/4;
    stb = pi/2;
    Va = 50; % 50 V on cone
    Vb = 0; % 0 V on plane ...it is ground
    
    % you want    ( 1/(sr^2 * sin(st) ) * d{st} [ sin(st) d{st} V ]  = 0
    % there is no change in sr or sf
    intg = c1/sin(st);
    V = int(intg, st);
    V = V + c2;
    c1s = Va / log(tan(sta/2));
    c2s = 0;
    V = subs(V, [c1, c2], [c1s, c2s]);
    V_pt = double(subs(V, st, pt_sph(2)));
    fprintf('at point, voltage is:  %.3f  V\n', V_pt);
    E = -1 * ee.getGradSph(V);
    E_pt = double(subs(E, [sr, st, sf], [pt_sph(1), pt_sph(2), pt_sph(3)]));
    fprintf('at point, E is:  [ %.3f, %.3f , %.3f ]  V/m\n', E_pt);
    % book is wrong
end


%------------------------------------------------------------------------------------------ #333
if sel == 333
    % solve Laplace's Equation
    temp1 = c1/sin(st);
    V_th1_th0 = int(temp1, st) + c2;
    
end

    
%{
    ee330 ch8 mag force, material, devices

    #1      ex1     physics
    #2      pp1     physics
    #3      ex2     physics
    #4      pe5     mag dipole
    #5      ex6     2 mag dipoles
    #6      ex7     magnetization
    #7      pp7     infer H,M,X_m given B
    #8      hw8,#2  fuck

%}
clc;
close all;
clearvars;
format short; % default, short, long, shortE, longE, shortG, longG, +, hex, rational [compact,loose]
format loose;


                sel = 8;  % CHANGE CHANGE CHANGE
                
                
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
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var

ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #1
if sel == 1
    m = 2; % mass in kg
    q = 3; % charge of C
    start = [1,-2,0];
    u = [4,0,3]; % m/s
    time = 1; % second
    E = [12, 10, 0]; % E field V/m
    
    % acceleration
    F = q .* E;
    a = double(F ./ m)  %  m/s^2
    
    % net velocity
    ux = subs(int(a(1), pt), pt, time);
    uy = subs(int(a(2), pt), pt, time);
    u_net = [ u(1) + ux , u(2) + uy, u(3) ]    %  [ 18t + 4, 15t , 3 ]
    
    % ke
    ke = (1/2) * m * norm(u_net)^2  % J
    
    % position
    pos = [ int(18*pt + 4, pt), int(15*pt, pt), int(3,pt) ];
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    m = 1; % kg
    q = 2; % C
    start = [0,0,0]; % and no intial velocity
    stop = [0,0,12];
    E = [0,0,3]; % V/m
    
    F = q .* E; % N
    a = F ./ m;
    u = [ 0 + 0, 0 + 0, 0 + 6];
    time = sqrt(4); % 2 s
    vf = [0,0,0 + time*a(3)]; % 12 m/s
    ke = (1/2)*m*vf(3)^2; % J
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    m = 2; % 2kg mass
    q = 1; % charge of 1 C
    u = [0,3,0]; % 3 m/s in y\hat
    B = [0,0,10];  % wb/m^2 mag field
    time = 4; % time to measure
    
    % calculate velocity and acceleration
    %   F= m d{u} = Q u x B
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    S = (1/10^4)*(10); % 10 cm^2
    I = 50; % A
    plane = 2*rx + 6*ry - 3*rz -7; % = 0
    g_plane = ee.getGradRec(plane);
    n_hat = g_plane ./ norm(g_plane);
    m = double((I*S) .* n_hat)
end


%------------------------------------------------------------------------------------------ #5
if sel == 5 % ex8.6
    m1 = [0,0,5]; % small current loop L1, has this momoent  Am^2
    p1 = [0,0,0]; % it is at the origin
    m2 = [0,3,0]; % another loop L2, has this moment Am^2
    p2 = [4,-3,10]; % located here
    % find torque on loop 2    T2 = cross(m2,B1)  ...due to B1
    m22 = ee.transVecRS(m2); % [3*sin(sf)*sin(st), 3*cos(st)*sin(sf), 3*cos(sf)]
    sin_st = 1/sqrt(5);
    cos_st = 2/sqrt(5);
    sin_sf = -3/5;
    cos_sf = 4/5;
    R = norm(p2);
    sub1 = [ sin(st), cos(st), sin(sf), cos(sf)];
    sub2 = [sin_st, cos_st, sin_sf, cos_sf];
    m222 = double(subs(m22, sub1, sub2))
end


%------------------------------------------------------------------------------------------ #6
if sel == 6 % ex8.7
    % inf slab, z[0,2]    
    ur = 2.5;
    B = (1e-3).*[10*ry, -5*rx,0]; % Wb/m^2
    J = double((1/(const.mu0*ur)) .* ee.getCurlRec(B)); %  -4,774.6  arz
    chim = ur-1;
    Jb = chim.*J; %  -7,162  arz
    M = (chim/(const.mu0*ur)) .* B
    % Kb to lower side of slab
end


%------------------------------------------------------------------------------------------ #7
if sel == 7 % pp8.7
    ur = 4.6;
    u = const.mu0*ur;
    chim = 1-ur; % 3.6
    B = [0,0,(1e-3)*10*exp(-ry)];
    H = B ./ u;
    Hval = double(subs(H,ry,0)); %  1730 exp(-ry)  z_hat
    M = (chim/u) .* B;
    Mval = double(subs(M,ry, 0))
end


%------------------------------------------------------------------------------------------ #8
if sel == 8 % hw8,#2
    syms rho_s0; assume(rho_s0, 'real');
    syms v0; assume(v0, 'real');
    syms a; assume(a, 'real');
    syms b; assume(b, 'real');
    syms x; assume(x, 'real');
    syms zp; assume(zp, 'real');
    syms z; assume(z, 'real');
    
    Rx = -x;
    Ry = (b/a)*sqrt(a^2 - x^2);
    Rz = (zp-z);
    R_yg0 = [Rx, -Ry, Rz];
    R_yl0 = [Rx, Ry, Rz];
    
    dS = sqrt( (b^2 - a^2)*x^2 + a^4 ) / ( a * sqrt(a^2 -x^2) ); % dx dz
    Rmag = ( sqrt( x^2 + (b^2/a^2)*(a^2 - x^2) + (zp-z)^2 ) )^3; 
    mult = dS / ( 4 * sym(pi) * Rmag );
    %pretty(mult);
    
    Jsx = ( rho_s0 * v0 * a * sqrt(a^2 - x^2) ) / sqrt( (b^2 - a^2)*x^2 + a^4 );
    Jsy = ( rho_s0 * v0 * b * x ) / sqrt( (b^2 - a^2)*x^2 + a^4 );
    Js_yg0 = [Jsx, -Jsy, 0];
    Js_yl0 = [Jsx, Jsy, 0];
    
    dH_yg0 = mult .* cross(Js_yg0, R_yg0);
    %pretty(dH_yg0);
    cp_yg0 = cross(Js_yg0, R_yg0);
    %pretty(simplify(cp_yg0));
    
    dH_yl0 = mult .* cross(Js_yl0, R_yl0);
    %pretty(dH_yl0);
    cp_yl0 = cross(Js_yl0, R_yl0);
    pretty(simplify(cp_yl0));
end

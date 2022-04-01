%{
     ee320 ch5 notes

    #1      pp5.1       transductance parameters
    #2      ex5.1       saturation w/ NMOS
    #3      pp5.2       min sat
    #4      ex5.2       nmos and regions
    #5      pp5.6       vds, id related in saturation
    #6      pp5.7       pmos
    #7      pp5.8       nmos design (find resistors)
    #8      pp5.9 /10      coupled NMOS
    #9      pp5.11      resistor is doubled to 13.1 k
    #10     pp14        pmos, resistor selection
    

    #5.561  hw10, 56a
    #5.562  hw10, 56b
    #5.582  hw10, 58b

    #991  lab8 

    

%}
clc;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');
format shortE; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]

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
global t; syms t; assume(t,'real'); % transform holder
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

                            sel = 5.582;  % CHANGE CHANGE CHANGE


%------------------------------------------------------------------------------------------ #1
if sel == 1
    L = .18e-6; % length
    tox = 4e-9; % thickness
    un = 450*(10^-4); %  m^2 / V s
    Vt = .5;
    Cox = const.ep_sox / tox;
    kn_p = un * Cox;                 %  388.13 uA/V^2
    
    % now find W so r_DS of 1k @ v_GS = 1
    r_DS = 1e3;
    v_GS = 1;
    v_OV = v_GS - Vt;
    temp = (un * Cox * v_OV)/L;
    temp1 = r_DS * temp;
    W = 1 / temp1;   % .93  um
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    L_min = .18e-6; % min length, for manufacturing process
    t_ox = 4e-9; % thickenss, 4 nm
    un = (1/10^4)*450; % m^2 / V s
    Vt = .5; % V
    
    C_ox = const.ep_sox / t_ox; %  8.6250e-03   F/m^2 
    kn_p = un * C_ox; %  3.8813e-04   A/V^2
    
    WL = 1.8/.18;
    kn = kn_p * WL; % 3.8813e-03
    id = 100e-6; % A   want to know vov, vgs, vdsm to operate in saturation
    
    vov = sqrt(2*id/kn); % 2.2700e-01 V
    vdsm = vov; % because in saturation...begins here
    vgs = Vt + vov; %  7.2700e-01 V
    
    % now find vov and vgs that make mosfet a 1k resistor for small vds
    rds = 1e3;
    vov = 1/(un*C_ox*WL*rds); %  2.5765e-01 V
    vgs = vov+Vt; % 7.5765e-01 V
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    vt = .35;
    wl = 10;
    id = 50e-6;
    L = 65e-9; % 65 nm process
    tox = 1.4e-9;
    cox = const.ep_sox/tox; %  2.4643e-02  F/m^2
    un = (1/10^4)*(216); 
    knp = un * cox; % 5.3229e-04
    kn = knp * wl;
    
    % in saturation:
    vov = sqrt((2*id)/kn); % 1.3707e-01  V
    vgs = vov + vt; % 4.8707e-01
    vdsm = vov; % min vds = vov if saturated  = .14 V
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    %tk = .18e-6; % thickness is .18 um
    L = .18; % in um
    W = 2; % in um
    Cox = 8.6e-15; % in fF/um^2
    un = (1/10^4)*450; % m^2/Vs
    Vtn = .5;
    
    kn_p = un * Cox*(10^12); % 387  uA/V^2
    kn = kn_p * (W/L); % 4.3 mA/V^2
    
    % find vgs, vds   on edge of saturation, id = 100 uA
    id = 100e-6;
    vov = sqrt((2*id)/kn); % .22 V
    vgs = Vtn + vov; % .72 V
    %vds = vov; % since on edge
    
    % vgs kept const .72V , id reduced ->  triode   id = 50 uA
    id = 50e-6;
    syms vds;
    func = (1/2)*vds^2 -vov*vds + id/kn;
    rotz = double(solve(func==0, vds))      % ignore vds > vov  choose vds = .06 V
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    kn_p = 400; % uA/m^2
    Va_p = 5; % V/um
    W = 16; % um
    L = .8; % um
    
    Va = Va_p * L; % 4 V
    lam = 1/Va; % .25  1/V
    
    vov = .2; % V
    vds = .8; % V
    kn = kn_p * (W/L);
    id = (1/2) * kn * vov^2 * (1 + vds*lam); % 192 uA
    id_p = (1/2) * kn * vov^2; % 160 uA
    ro = Va / (id_p/10^6); % 25 k ohm
    
    % now vds increases by 1 V
    vds = vds + 1; % 1.8 V
    id1 = (1/2) * kn * vov^2 * (1 + vds*lam) % 232 uA
    chng = id1-id; % up 40 uA
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    Vtp = -.5; % V
    kp_p = 100; % uA/V^2
    WL = 10; % W/L
    vs = 1.8;
    kp = kp_p * WL / 10^6;
    
    % find range vg PMOS conducts ...must be lower than than source by abs(Vtp)
    vg = vs-abs(Vtp); % anything less than 1.3 V and it will turn on/conduct   vg < 1.3V
    % vd >= vg + abs(Vtp)  and it is in triode                                 vd >= vg + .5
    % vd <= vg + abs(Vtp) and saturated                                        vd <= vg + .5
    
    %neglect modulation (lam=0) find abs(vov) and vg corresponding to range vd , id = 50 uA
    id = 50e-6;
    vov_temp = sqrt((id*2)/kp); % .32 V
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    Vtn = .4; % V
    kn_p = 400/10^6; % A/V^2
    WL = 5/.4; % W/L
    kn = kn_p * WL;
    Id = 100e-6;
    Vd = .2; % V
    Vdd = 1; % V
    Vss = -1; % V
    Vg = 0;
    
    Rd = (Vdd-Vd)/Id; % 8 k oms
    vov = sqrt((Id*2)/kn);
    Vs = Vg - Vtn - vov;
    Rs = (Vs-Vss)/Id; % 4 k ohms
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    Vdd = 1.8; % V
    Vd = .7; % V
    Vtn = .5; % V
    kn_p = .4/10^3; % A/V^2
    WL = .72/.18; % W\L
    kn = kn_p * WL;
    
    R = 2*(Vdd-Vd)/(kn*(Vd-Vtn)^2); % 34.4 k ohm  had to infer
    id = (Vdd-Vd)/R; % 32 uA
    
    Vd2 = Vd - Vtn; % at edge of saturation
    Vgs2 = Vd; % since connected and sources are grounded
    id2 = (1/2)*kn*(Vgs2-Vtn)^2;
    R2 = (Vdd-Vd2)/id2; % 50 k ohm
    
    %R2 = 2*(Vdd-Vd2)/(kn*(Vgs2-Vtn)^2)
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    kn = .002;
    Rd = 13.1e3; 
    Vdd = 2;
    Vg = Vdd;
    Vs = 0;
    Vtn = .5;
    % must be triode, vgs > vtn
    Vov = Vg-Vs-Vtn;
    syms x;
    polyn = (-1/4)*Rd*kn*x^2 + x*(1+(1/2)*Rd*kn*Vov) - Vdd;
    rotz = double(solve(polyn==0, x)); % use Vd = .1
    Vd = rotz(1,1)
    Id = (Vdd-Vd)/Rd; % about .15 mA 
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    Vss = 1.8;
    vov = -.6;
    vtp = -.4;
    kp = (.1e-3)*(10/.18);
    % they imply saturation.... Vdg<abs(Vtp)  
    
    vs = 1;  %  vsg = abs(vtp)+abs(vov)    vs-vg = .4 + .6 = 1
    id = (1/2)*kp*vov^2;
    R = (Vss-vs)/id;
end
    

%------------------------------------------------------------------------------------------ #5.561
if sel == 5.561
    syms v2;
    lhs = (v2+1)/4000;
    rhs = (1/2)*(.005)*(-v2-.4)^2;
    %eqn = rhs-lhs;
    %pretty(simplify(eqn))
    rotz = double(solve(lhs==rhs, v2));
    v2 = rotz(1,1);
    
    syms v1;
    v1eqn = (1/2)*(.005)*(.1-v1)^2 - .0001;
    rotz = double(solve(v1eqn==0, v1))
    v1 = rotz(1,1);
end


%------------------------------------------------------------------------------------------ #5.562
if sel == 5.562
    syms v5;
    lhs = v5/4000;
    rhs = (1/2)*(.005)*(.6-v5)^2;
    rotz = double(solve(lhs==rhs,v5))
end


%------------------------------------------------------------------------------------------ #5.582
if sel == 5.582
    syms v4;
    lhs = (1/2)*(.002)*(v4-.4)^2;
    rhs = (1/2)*(.0005)*(1.6-v4)^2;
    rotz = double(solve(lhs==rhs,v4))
end
    
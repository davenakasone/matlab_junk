%{
    ee320 ch4

    #1  4.3 pp  basic diode models
    #2  4.1 ex
    #3  4.4 ex
    #4  4.3 ex   p216  forward iv 
    #5  4.7 ex   p218
    #6  4.8 ex   p218
    #7  4.9 ex   p218
    #8  4.10 ex  p218
    #9  4.11 ex  p219   reverse iv
    #10 4.4 ex   p221  iterative approch to exp model
    #11 4.16 pp  p226   zener diode
    #12 4.5 ex   p230   small voltage estmation
    #13 4.17 ex  p233    ""
    #14 4.18 ex  p233
    #15 4.19 pp  p234 design with small signal
    #16 4.20 ex  p236   zener on small signal
    
    
    
%}
clc;
close all;
clearvars;
sympref('MatrixWithSquareBrackets', 1);
sympref('PolynomialDisplayStyle', 'ascend');
    old_val = sympref('HeavisideAtOrigin', 1);


                sel = 921;  % CHANGE CHANGE CHANGE


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
global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global ik; syms ik; assume(ik, {'real', 'integer'}); % index k, usually for convolution
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
global T; syms T; assume(T, 'real');  % period bound, continous case
global t0; syms t0; assume(t0, 'real'); % arbitrary t0, usually offset
global tau; syms tau; assume(tau, 'real'); % dummy tau, usually for convolution

ee = cls_EE330_helper();
const = cls_CONST();
%publisher('ee330_hw4.m');
%------------------------------------------------------------------------------------------ #1
if sel == 1
    vi = 10; %pp
    R = 5e3; % ohm
    pd = 2*pi/2; % a period of normally 2 pi is halfed
    id = vi / R % 0020  A
    dc = vi / pd;  % 3.1831 V
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    vs = 24; % 24v peak amp
    bat = 12; 
    R = 100;
    % it has to exceed battery's 12V to conduct
    % 24(cos[th]) = 12     cos([th]) = 1/2   th = pi/3
    th = pi/3;
    cond_angle = 2*th; % 2.0944   1/3 of a cycle it is actually conducting
    id = (vs - bat)/R;  % .1200 A
    max_rev = vs + bat; % account for negative source and reverse 36 = 24 + 12
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    i_a = 5/(2.5e3); % 0020 A
    v_a = 0; % voltage already expended
    i_b = 0;   % revese blocks
    v_b = 5;   % no change
    i_c = 0; % blocked, needs positive polarity
    v_c = 5; % no change   remember 0- -5 = 5
    i_d = -5/(2.5e3); % 0.0020 A  going through diode...reverse arrow
    v_d = 0; % potential expended
    i_e = (3)/(1e3); % you take the largest input..it will negate others, reversing/stopping them
    v_e = (3); % only take the source that doesnt reverse the rest
    i_f = (1-5)/(1e3); %0.0040  A  this stops other 3
    v_f = 1; % only voltage making it...or ohm's law remainder
    
    i_avg = .001; % half cycles are pi  not 2 pi
    r_meter = 50;
    v_in = 10; % 20 peak to peak
    maxr = ((v_in/pi)/i_avg)-50; % 3.1331e+03  ohms
    
    r20 = (20-5)/.05;  %max 20V  just subtract internal source and ohm   r=300
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    % 1 mA rated diode displays fwd V = .7    if i = 1 mA
    % 1 A diode made to conduct at .7 range
    id_1ma = .001; % diode current for the 1 mA rated diode
    vf = .7; % designed to exibit behvior in this range
    id_1a = 1; % diode current for 1 A, higher powered diode
    Is = id_1ma * exp( -vf / const.sil_vtRoom); % Is must be  6.9144e-16  A
    % the little diodes are a 1:1 scaling const of big guy, so 1000x for Is
    Is1000 = 1000*Is;  % 6.9144e-13  A
    % or
    is = id_1a * exp(-vf/const.sil_vtRoom); % 6.9144e-13
end

%------------------------------------------------------------------------------------------ #5
if sel == 5
    i_beg = .1e-3; % current is .1 mA to begin
    i_end = 10e-3; % current is 10mA after some temperature change
    dV = const.sil_vtRoom * log(i_end/i_beg); % 0.1151 V
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    v = .7;
    i = 1e-3;
    i1 = .1e-3;
    i2 = 10e-3;
    Is = i * exp(-v/const.sil_vtRoom); % find this first, and everything else is good
    
    v1 = const.sil_vtRoom * log(i1/Is);  %  0.6424 V
    v2 = const.sil_vtRoom * log(i2/Is); % 0.7576 V
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    v = .3; % a germanium diode's voltage drop
    i = .2e-3;
    Is = i * exp( -v / const.sil_vtRoom); % not sure about VT, but 1.2288e-09
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    Is = 10^-14;
    kel = const.cel2kel(25);
    Is_125 = Is*(1.15)^(125-25); %  1.1743e-08 A    increase 15% every 1C° every
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % high current diode where revese leakage independt of voltage
    % V = 1 @ 20°C, find V at 40°C and 0°C
    vs = 9;
    v = 1;
        %Vt20 = const.sil_vtRoom;   not needed
        %vt40 = const.cel2vt(40);
        %vt0 = const.cel2vt(0);
    R = 1e6;
    i = v/R;
    Is = -i;  % based on relationship 20°C
    
    Is40 =(Is*2)*2 % doubles every 10°C
    i40 = -1*Is40;
    V40 = Is40 * R; % keep direction, 4 V
    
    Is0 = (Is/2)/2; % halves every 10°C
    i0 = -1*Is0;
    V0 = i0 * R;  
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    % want ID and VD, VDD = 5 V, R = 1k, assume 1mA @ .7V
    VDD = 5;
    R = 1e3;
    i = 1e-3;
    VDa = .7; % begin by assuming VD = .7 V
    IDa = (VDD - VDa)/R; % implies 0.0043 A
    % use that log eqn
    V2a = VDa + (2.3 * const.sil_vtRoom * log10(IDa/i)); % 0.7364
    
    % iterate again
    IDb = (VDD - V2a)/R; % 0.0043 A
    V2b = V2a + (2.3 * const.sil_vtRoom * log10(IDb/IDa)); % 0.7362   close enough, no more iteration
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    safe_p = .2; % W
    Vz = 3.5; % V
    It = .01; % test current of 10 mA
    r_z = 10; % incremental resistance
    
    i = .02; % 20 mA through
    v20 = Vz - (It-i)*r_z; %  3.6000 V
    r = v20/i; % ohm  180
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    vs = 10; % also has ac comp cos(377t)
    vs_ac = 1;% % ac comp has 1 v max amplitude
    R = 10e3;
    vd = .7; % assuming a .7 V drop
    id = 1e-3; % this drop occurs at 1 mA
    
    ID = ( vs - vd ) / R; % calculate DC curret through diode = 9.3000e-04 A
    % this is close to the given 1mA, so voltage will be close to assumed .7
    
    rd = const.sil_vtRoom / ID; % linearize circuit with this resistor = 26.8817 ohms
    
    vo = vs_ac * ( rd / (R + rd) ); % find ac comp = 0.0027 V    it's good, < 5 mV  
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    % small signal resistance  ID = .1 mA, 1 mA, 10 mA
    rd1 = const.sil_vtRoom / .1e-3; % 250
    rd2 = const.sil_vtRoom / 1e-3; % 25
    rd3 = const.sil_vtRoom / 10e-3; % 2.5
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    dv = -10e-3; % given change in voltage
    ID = 1e-3; % given bias = 1 mA
    rd = const.sil_vtRoom / ID; % find equiv resistance
    did = dv / rd; % implied change in current    -4.0000e-04 A
    expID = ID * ( exp( dv / const.sil_vtRoom ) - 1 ); % -3.2968e-04 A   close
    
    dv = -5e-3; % given change in voltage
    ID = 1e-3; % given bias = 1 mA
    rd = const.sil_vtRoom / ID; % find equiv resistance
    did = dv / rd; % implied change in current    -2.0000e-04 A
    expID = ID * ( exp( dv / const.sil_vtRoom ) - 1 ); % -1.8127e-04 A   close
    
    dv = 5e-3; % given change in voltage
    ID = 1e-3; % given bias = 1 mA
    rd = const.sil_vtRoom / ID; % find equiv resistance
    did = dv / rd; % implied change in current    2.0000e-04 A
    expID = ID * ( exp( dv / const.sil_vtRoom ) - 1 ); % 2.2140e-04 A   close
    
    dv = 10e-3; % given change in voltage
    ID = 1e-3; % given bias = 1 mA
    rd = const.sil_vtRoom / ID; % find equiv resistance
    did = dv / rd; % implied change in current    4.0000e-04 A
    expID = ID * ( exp( dv / const.sil_vtRoom ) - 1 ); % 4.9182e-04 A   close
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    % want vo = 3 when IL = 0, vo to change 20mV per 1 mA of load current
    vs = 15;
    dv = 20e-3;
    di = 1e-3;
    r = dv/di; % implies 20 ohms for all 4 diodes  (vo spans the stack)
    rd = r/4; % 5 ohms per diode
    
    ID = const.sil_vtRoom / rd;  % each diode has ID = 0.0050 A
    R = (vs - 3) / ID;% 2400  by voltage division
    Is = ID * exp( -.75 / const.sil_vtRoom ); %  4.6788e-16  A   each loses 3/4 volts
    VD = const.sil_vtRoom * log(4*di/Is) % 0.7444 V  per diode
    VD1 = VD*4; % 2.9777 V out
    dif =  VD1 - 3;  % decrease  -0.0223 V
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    vz = 6;
    iz = 5e-3;
    rz = 80;
    
    % iz halfed means voltage has to reduce
    di = (iz/2)-iz;
    dv = di*rz; % (decrease) -0.2000 V
    v = vz + dv; % 5.8000
    
    % iz doubled means voltage has to increase
    di = (2*iz)-iz;
    dv=di*rz; % (increase) 0.4000 V
    v = vz + dv; %  6.4000
    
    vz0 = vz - (iz*rz); % 5.6000 V
end

%------------------------------------------------------------------------------------------ #17
if sel == 17
end




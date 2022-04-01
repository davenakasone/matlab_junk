%{
    ch5, E in material space    notes

    #1      ex5.1       given J, find I    sph
    #2      pp5.1       given J, find I    cyn
    #3      ex5.2       van de graaff generator
    #4      pp5.2       " "
    #5      ex5.3       wire and electrons
    #6      pp5.3       use J,I,sig to derive others
    #7      ex5.4       finding resistence
    #8      pp5.4       this time 2 materials, find resistance ...parallel resistors
    #9      ex5.5       dielectric cube, radial polarization, find bound charge densities --> 0
    #10     pp5.5       long rod, find bound charges and show net 0
    #11     ex5.6       polystyrene and some basic E, D, P, ect
    #12     pp5.6       paralell plate capacitor...more basics
    #13     ex5.7       sphere and charge
    #14     pp5.7       chi, D, E, given  / unknowns
    #15     ex5.9       dielectric:dielectric interface  apply boundary conditions
    #16     pp5.9       dialectric : free space
    #17     ex5.10      conductor : dialectric
    #18     pp5.10      conductor : free space



%}
clc;
close all;
clearvars;


                sel = 18;  % CHANGE CHANGE CHANGE


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


ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #1
if sel == 1
    J = (1/sr^3) .* [ 2*cos(st), sin(st), 0];
    %  a)  find I, sr = .2, st[0,pi/2] , sf[0,2 pi]    I = int(int( dot(J, dS), st, sf...
    % <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    >
    dS = [ sin(st)*sr^2, 0, 0];
    intg = dot(J, dS);
    intg = subs(intg, sr, .2);
    I = int(int(intg, st, 0, sym(pi)/2), sf, 0, 2*sym(pi));
    fprintf('a) the current is:  %.3f  A\n', I);
    %  b)  I, shell sr = .1
    intg = dot(J, dS);
    intg = subs(intg, sr, .1);
    I = int(int(intg, st, 0, sym(pi)), sf, 0, 2*sym(pi));
    fprintf('b) the current is:  %.3f  A\n', I);
    chk = ee.getDivgSph(J); % divg(J) = 0
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    J = [10*cz*sin(cf)^2, 0, 0];
    % <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} > 
    dS = [ cr, 0, 0];
    intg = dot(J, dS);
    intg = subs(intg, cr, 2);
    I = int(int(intg, cf, 0, 2*sym(pi)), cz, 1, 5);
    fprintf(' I = %.3f  A\n', I);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    pS = (1e-7); % C/m^2
    u = 2; % m/s
    t = 5; % seconds
    width = .1; % m
    I = pS * u * width; % current on dome is charge/surf area  * surface area * velocity
    Q5 = I * t; % charge collected over 5 sec
    fprintf('over 5 seconds, charge is: %.3f  nC\n', Q5*(1e9));
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    width = .1;
    u = 10;
    ohms = 1e14;
    pS = .5e-6;
    % when steady, current through leakage path = charge transposted / unit time
    I = pS * u * width; % .. V = IR
    V = I * ohms;
    fprintf('potential difference is %.3f  MV\n', V*(1e-6));
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    diameter = .001; % 1mm
    sig = 5*10^7; % S/m
    ne = 10^29; % free electrons / m^3
    E = .01;  % apply field of 10mV
    % a) charge density of free electrons  pV
    pV = ne * -const.ef;
    fprintf('charge density of free electrons is pV = %.3f 10^10 C/m^3\n', pV*(10^-10));
    % b) current density  J = sig E
    J = sig * E;
    fprintf('current density J = %.3f  kA/m^2\n', J*(1e-3));
    % c) wire current I = JS
    I = J * pi * ( diameter/2 )^2;
    fprintf('current I in wire = %.3f\n', I);
    % d) drift velocity of electrons
    u = J/pV;
    fprintf('drift velcotity of electrons is %d  m/s\n',u);
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    pV = 1.81e10; % free charge density in copper  C/m^3
    J = 8e6; % A/m^2
    sig = const.cond_cop; % conductivity of copper, S/m
    n_free = pV/const.e; % find free electrons per unit volume
    u = J/pV; 
    E = J/sig;
    fprintf(' drift velocity is %d  m/s   E = %.3f V/m\n', u, E);
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    sig = const.cond_led; % S/m  given conductivity
    leng = 4; % 4m in length
    % the cross section S is uniform  R = l/ (sig * S)
    % S = d^2 - pi(r^2)
    r = .01/2;
    S = .03^2 - pi*r^2;
    R = leng / (sig * S);
    fprintf('resistance is %.2f micro Ohms\n', R*(1e6));
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    sig_ld = const.cond_led; % S/m for lead
    sig_co = const.cond_cop; % S/m for copper
    leng = 4; % 4m in length
    
    r = .01/2;
    d_ld = .03; 
    S_ld = d_ld^2 - pi * r^2; % cross section of lead material
    S_co = pi * r^2; % cross section of copper material
    
    R_ld = leng / ( sig_ld * S_ld );
    R_co = leng / ( sig_co * S_co );
    R = const.para(R_ld, R_co);
    fprintf('total resistance, end to end, is %.2f micro Ohms\n', R*(1e6));
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    syms a; assume(a,'real');
    syms L; assume(L, 'real'); assume(L>0);
    P = a .* [rx, ry, rz];  % P = ar   radially outward
    
    % for each of the 6 faces, there is a charge density rho_PS
    % the cube sits at orgin, x[L/2 , -L/2] , y[L/2 , -L/2], z[L/2 , -L/2]
    rho_ps1 = dot(P, [1, 0, 0]);
    rho_ps1 = subs(rho_ps1, rx, L/2); % (L*a)/2
    rho_ps2 = dot(P, [-1,0,0]);
    rho_ps2 = subs(rho_ps2, rx, -L/2); % (L*a)/2
    rho_ps3 = dot(P, [0,1,0]);
    rho_ps3 = subs(rho_ps3, ry, L/2); % (L*a)/2
    rho_ps4 = dot(P, [0,-1,0]);
    rho_ps4 = subs(rho_ps4, ry, -L/2); % (L*a)/2
    rho_ps5 = dot(P, [0,0,1]);
    rho_ps5 = subs(rho_ps5, rz, L/2); % (L*a)/2
    rho_ps6 = dot(P, [0,0,-1]);
    rho_ps6 = subs(rho_ps6, rz, -L/2); % (L*a)/2
    
    % total bound surface charge is rho_ps over the faces, they are all the same, so that is good
    % just multiply by 6
    Qsb = int(int(6 * rho_ps1, ry, (-L/2), (L/2)), rz, (-L/2), (L/2)); % 3*a*L^3
    
    % volume charge density is -divg(P)
    rho_pv = -1 * ee.getDivgRec(P); % -1 * ( a + a + a ) = -3*a
    Qvb = int(int(int(rho_pv, rx, (-L/2), (L/2)), ry, (-L/2), (L/2)), rz, (-L/2), (L/2)); % -3*a*L^3
    
    Qt = Qsb + Qvb; % = 0 as expected
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    syms L; assume(L, 'real'); assume(L>0); % rod is on x axis, [0,L]
    syms a; assume(a, 'real');
    syms b; assume(b, 'real');
    syms A; assume(A, 'real'); assume(A>0); % cross-section area
    P = [ b + a * rx^2, 0, 0 ]; % given polarization
    
    rho_ps1 = dot(P, [-1,0,0]);
    rho_ps1 = subs(rho_ps1, rx, 0); % -b
    rho_ps2 = dot(P,[1,0,0]);
    rho_ps2 = subs(rho_ps2, rx, L); % b + a * L^2
    % could integrate, or just multiply by cross sectional area
    Qsb = rho_ps1 * A + rho_ps2 * A; % A*(b + L^2*a) - A*b    =  a*A*L^2

    % for volume charge density
    rho_pv = -1 * ee.getDivgRec(P); % -2*a*rx
    rho_pv0 = subs(rho_pv, rx, 0); % 0  at x=0
    rho_pvL = subs(rho_pv, rx, L); % -2*L*a   at x = L
    % be smart and integrate area over one integral
    Qvb = int( A * rho_pv, rx, 0, L); % -a*A*L^2
    
    Qt = simplify(Qsb + Qvb); % 0 as expected
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    epr = const.epr_pst; % given polystyrene
    Emag = 10e3; % strength is given as 10 kV/m
    d = 1.5e-3; % parallel plates have distance of 1.5mm
    
    Dmag = double(Emag * const.ep0 * epr)*(1e9); %  225.4695  nC / m^2
    Pmag = double(Emag * (epr-1) * const.ep0)*(1e9); % 137.0501 nC / m^2
    rho_S = Dmag; % because rho_S (free) = dot(D,n)    +/-  on each plate  all same direction
    rho_ps = Pmag; % because rho_ps = dot(P,n)   +/- on each plate, same direction
    V = Emag*d; % 15 V   also, E = V/d  
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    d = 2e-3; % plates separated by 2mm     assume x = 0, x = .002 for plate location
    V = 1e3; % 1 kV is applied
    epr = const.epr_pst; % polystrene is between plates = 2.55
    
    E = [V/d, 0, 0]; % 500 kV/m from x = 0 to x = .002
    temp = (epr-1)*const.ep0;
    P = double(temp .* E); %  [6.853 ,0,0] uC/m^2
    
    rho_ps1 = dot(P, [1,0,0]);
    rho_ps1 = double(subs(rho_ps1, rx, d)); % 6.8525  uC/m^2
    rho_ps2 = dot(P,[-1,0,0]);
    rho_ps2 = subs(rho_ps2, rx, 0);
    rho_ps = rho_ps1 + rho_ps2
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    epr = 5.7; % given dielectric const
    r = 10e-2; % 10cm radius
    q = 2e-12; % 2 pC put in center of sphere
    qex = -4e-12; % a -4 pC charge is put on the surface of the sphere
    
    % rho_ps surface density of polarization charge on sphere's surface
    % assume charge is at origin and use Gauss or Coulumb
    E = [  q / (4 * pi * epr * const.ep0 * sr^2) , 0 , 0 ];
    E = double(subs(E, sr, r)); % [ .32,0,0] V/m
    P = double(((epr-1) * const.ep0) .* E); % [13.12,0,0] pC/m
    rho_ps = dot(P,[1,0,0]); % 13.12 pC/m^2
    
    %by coulumb's law, sphere is positive, so it will try to pull in external charge
    F = [double(( q * qex ) / ( 4 * pi * const.ep0 * epr * r^2 )) , 0 , 0 ]; % [ -1.26,0,0] pN
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    Ex = 5; % V/m  missing Ey, Ez
    P = ( (1e-9) * (1/(10*pi)) ) .* [3, -1, 4] ; % nC/m^2
    
    chi = double(P(1) / ( const.ep0 * Ex ));  % 2.1600
    epr = 1 + chi; % 3.1600
    ep = double(const.ep0 * epr); % 2.7941e-11
    Ey = double(P(2) / ( const.ep0 * chi ));
    Ez = double(P(3) / (const.ep0 * chi ));
    E = [ Ex, Ey, Ez ];
    D = (ep*(1e12)) .* E; % [ 139.7027  -46.5676  186.2702 ]  pC/m^2
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    % 2 homo, iso, linear  (simple) dielectrics meet on z=0 plane
    ep_r1 = 4; % z>0 
    ep_r2 = 3; % z<0
    E1 = [5e3, -2e3, 3e3]; % kV/m  z >= 0
    n1 = [0,0,1];  % arz is normal to boundary plane
    E1n_mag = dot(E1, n1);
    E1n = [0,0,E1n_mag];
    E1t = E1 - E1n; % implied tangent
    E2t = E1t; % because of boundary condition
    D1n = ep_r1 * E1n;
    D2n = D1n; % second boundary condition
    E2n = (ep_r1/ep_r2).*E1n;  % or ep_r2 * D2n
    
    E2 = E2t + E2n; % easy after set up
    tan_th1 = norm(E1t) / norm(E1n);  % or  dot(E1,n) = norm(E1) * 1 * cos(th1) and solve
    th1 = rad2deg( atan( tan_th1) );
    tan_th2 = norm(E2t) / norm(E2n); % or solve by dp
    th2 = rad2deg(atan( tan_th2));
    alp1 = 90-th1;
    alp2 = 90-th2; % the alpha angles are what you want...they determine interface
    chk1 = (tan(deg2rad(th1))/tan(deg2rad(th2))) - (ep_r1/ep_r2);
    
    %enegery desity is  we = 1/2 ep norm(E)^2
    we1 = double((1/2)*(const.ep0 * ep_r1)*norm(E1)^2)*(1e6);  % uJ/m^3
    we2 = double((1/2)*(const.ep0 * ep_r2)*norm(E2)^2)*(1e6);  % uJ/m^3
    
    % if there is a cube with sides = 2m, centered at (3,4,-5)
    % then x[2,4] y[3,5] z[-6,-4]  it is entirely in second region
    energy = we2 * 2 * 2 * 2; % because it is a energy volume density * actual volume
    energ = double(int(int(int(we2, rx, 2,4),ry,3,5),rz,-6,-4)); % or integrate
    chk2 = energy - energ
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    epr1 = 2.5;
    D1 = [ 12, -10, 4]; % nC/m^2
    D1n = [ D1(1), 0, 0]; % what is normal to boundary
    D1t = D1 - D1n; % best, just to be sure
        %D1t = [ 0, D1(2), D1(3)]
   D2n = D1n;
   
   % Dt2/1 = Dt1/epr1    Dt2 = 1/epr1 Dt1
   Dt2 = (1/epr1) .* D1t;
   Dt2 = Dt2 + D2n; % nC / m^2
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    % y < 0 is a perfect conductor
    % y > 0 is a dialectric e1r = 2
    q = 2e-9; % 2nC surface charge on conductor
    epr = 2; % for conductor...free space, epr = 1, ep = ep0
   
    % @ (3, -2, 2) you are inside conductor, E=D=0
    % @ (-4, 1, 5) you are in dielectric    Dn = rho_s = q (given)
    Dn = [ 0, q, 0]; % normal from conductor to dielectric    must be completely normal to surface
    D = Dn; % conductor, by properties
    E = double(D ./ (const.ep0 * epr )); % v/m
    
    % just use and abuse conductor properties
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    E = (1e-3).*[60, 20, -30]; % mV/m  at a particular point on interface of air and conductor
    % this must be completley normal on the bondary, by properties of conductor
    D = double( (const.ep0*(1e12)) .* E); %pC/m^2
    rho_s = norm(D); % pC/m^2
end
    
    
    
    
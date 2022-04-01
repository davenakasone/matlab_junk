%{
    chapter 4 problems

    #1      basic coulomb, 2 point charges
    #2      manipulation of coulomb
    #4      equate electrical force to levitate
    #5      determine total charge
    #7      volume density --> charge in cylinder wedge
    #8      surface density over triagular region   ...split it up
    #9      wedge shaped surface...rec triangle
    #10     wedge volume cyn, find charge
    #11     rec line and finding E it produces
    #12     ring/loop exerts a force
    #14     infinite sheet
    #15     plane not on axis carrying charge
    #17     3 surfaces...closest point, super position for E
    #21     ring charge, find D
    #22     maxwell 1  D --> rho_v
    #23     flux from implied E   .... D = ep0 E 
    #24     D --> rho_V by maxwell   not sure which one to use
    #25     charge desnity pV from D
    #26     str8 Gauss
    #27     2 spheres, use gauss
    #29     maxwell 1  --> Guass    rec
    #30     "                 "     cyn
    #33     net flux on sph, but watch the limits of pV ...forced to work surface
    #35     potential , 2 pnt charges
    #36     potential , 4 point charges   poi above
    #37     potential as charge spaces around loop
    #38     potential with a few point charges on a poi
    #39     static Maxwell 2   V = -E
    #40     V --> E (maxwell_1) --> Qenc by maxwell_2
    #43     given D, find pV @ pnt
    #44     shells of sphere
    #45     pV from D, maxwell 1
    #46     work to move charge, given 2pts  rec, 3 different ways
    #47     work cyn
    #48     work sph  point charges
    #49     potential for points in sph
    #50     had to see what way to convert
    #51     sheet and work
    #52     potential --> E --> D -->pV
    #53     potential --> E, prove conservative
    #54     sph, V-->D
    #57     genuine flux
    #60     dipole
    #62     dipole
    #64     energy/work for charge transfer
    #65     volume enegery box
    #66     volume enegery hemisphere
    #68     E , find energy in volume
    #69     V --> E, energy in volume


%}
clc;
close all;
clearvars;


                sel = 69;  % CHANGE CHANGE CHANGE
                

global rx; global ry; global rz; % rectangular params  rx, ry, rz        
global cr; global cf; global cz; % cylindrical params  cr, cf, cz
global sr; global st; global sf; % spherical params    sr, st, sf
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); syms cf; assume(cf,'real'); syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); syms st; assume(st,'real'); syms sf; assume(sf, 'real');
global pL; syms pL; assume(pL, 'real'); % line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); % surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); % volume charge denisty in C / m^3
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve


ee = cls_EE330_helper();
const = cls_CONST();
%------------------------------------------------------------------------------------------ #1
if sel == 1                             % [ -1.8291 , -0.5226 , 1.3065 ]  mN
    q1 = 5e-6;
    r1 = [3, 2, 1];
    q2 = -4e-6;
    r2 = [-4, 0, 6];
    
    % want " force on q1   F21 --> " force charge 2 puts on charge 1 "
    R = r1 - r2; % poi - source  =  direction from q2 to q1
    F21 = (q1 * q2 / ( 4 * sym(pi) * const.ep0 * norm(R)^3 ) ) .* R;
    fprintf(' force on q1 due to d2 is F = [ %.4f , %.4f , %.4f ]  mN\n', F21*(1e3));
    
    F = ee.CLFpt(q1, r1, q2, r2); % obj method works
    fprintf(' force on q1 due to d2 is F = [ %.4f , %.4f , %.4f ]  mN\n', F*(1e3));
end


%------------------------------------------------------------------------------------------ #2
if sel == 2                            % q1 = -8.3232  nC   ,  q1 = -44.9452  nC
    q1 = 0; % need to find
    r1 = [4, 0, -3];
    q2 = 4e-9;
    r2 = [2, 0, 1];
    
    % want E with arz = 0
    rE = [5, 0, 6];
    r2_E = rE - r2; % q2 to poi direction vector
    E2 = (q2 / ( 4 * sym(pi) * const.ep0 * norm(r2_E)^3) ) .* r2_E; % Arz must cancel with q1
    fprintf('have to make z componenet 0  --> E2 = [ %.4f , %.4f , %.4f ]  V/m\n', E2);
    fprintf('q1 must produce E1 = [ ? , ? , 0 ]  V/m\n');
    
    r1_E = rE - r1; % q1 to poi direction vector
    q1 = -q2 * (5/34^1.5) * (82^1.5/9);
    E1 = (q1 / ( 4 * sym(pi) * const.ep0 * norm(r1_E)^3) ) .* r1_E;
    E = double(E1 + E2);
    fprintf(' with q1 = %.4f  nC , E = [ %.4f , %.4f , %.4f ]  V/m\n', (1e9)*q1, E);
    
    % now a test charge at same point will have x = 0 component
    tst = [ 5, 0, 6];
    Qtst = 1; % for convienence won't matter, cancels
    q2_tst = tst - r2; % direction vector of q2 to the test charge
    F_q2_tst = (q2 * Qtst / ( 4 * sym(pi) * const.ep0 * norm(q2_tst)^3)) .* q2_tst;
    
    q1_tst = tst - r1; % q1 to test charge
    q1 = -q2*(3/34^1.5)*(82^1.5);
    F_q1_tst = ( q1 * Qtst / ( 4 * sym(pi) * const.ep0 * norm(q1_tst)^3)) .* q1_tst;
    F = double( F_q1_tst + F_q2_tst);
    fprintf('\n q1 = %.4f  nC   produce F = [ %.4f , % .4f , %.4f ]  N\n',q1*(1e9), F); 
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    m = 2; %kg
    g = 9.8; % m/s^2
    f = m*g;
    q = -4e-3; % -4 mC
    E = f / q; % because F = q E
    fprintf('E must be [ 0, 0, %.4f ]  V/m   \n', E);  % [ 0, 0, -4900.0000 ]  V/m  
end


%------------------------------------------------------------------------------------------ #5
if sel == 5                         %  0.500  C      1.2064  mC         1579.1367  C
    
    % line x[0,5]
    pL = [ (1e-3)*12*rx^2 , 0, 0]; % given density
    dl = [1, 0, 0];
    intg = dot(pL, dl);
    Q = int(intg, rx, 0, 5);
    fprintf(' this line charge has a total charge of:  %.3f  C\n', Q);
    
    % cylinder  cr = 3 , cf[0, 2 pi] , cz[0,4]
    pS = (1e-9) * cr * cz^2; % given
    % sides  cr = 3 , cf[0, 2 pi] , cz[0,4]     cr d{cf} d{cz}
    intg = pS * cr;
    intg = subs(intg, cr, 3);
    Q_sides = int(int(intg, cr, 0, 2*sym(pi)), cz, 0, 4);
    % top   cr[0,3] , cf[0, 2pi] , cz = 4      cr d{cr} d{cf}
    intg = pS * cr;
    intg = subs(intg, cz, 4);
    Q_top = int(int(intg, cr, 0, 3), cf, 0, 2*sym(pi));
    % bottom   cr[0,3] , cf[0, 2pi] , cz = 0      -cr d{cr} d{cf}
    intg = pS * -cr; % orientate
    intg = subs(intg, cz, 0);
    Q_bot = int(int(intg, cr, 0, 3), cf, 0, 2*sym(pi));
    Q = Q_sides + Q_top + Q_bot;
    fprintf(' this cylinder surface has a charge of %.4f  mC\n', Q_sides*(1e6)); % just wants sides
    
    % sphere, sr = 4 , st[0,pi] , sf[0, 2pi]
    pV = 10 / (sr * sin(st));
    intg = pV * sin(st) * sr^2;
    Q = int(int(int(intg, sr, 0, 4), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    fprintf(' this sphere volume has a charge of %.4f  C\n', Q);
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    % volume is cr[0,2] , cf[pi/6 , pi/2] , cz[0,1]  a cyn wedge
    pV = (1e-3)*5*cz*cr^2; % given vol density C/m^3
    intg = pV * cr; % pV * cr d{cr} d{cf} d{cz}
    Q = int(int(int(intg, cr, 0, 2), cf, sym(pi)/6, sym(pi)/2), cz, 0, 1);
    fprintf('this wedge has a total charge of: %.4f  mC\n', 1000*double(Q)); % 10.4720  mC
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % rz = 0
    % part 1:  rx[0,2] , ry[0, rx]    y=0 to y = x
    % part 2:  rx[2,4] , ry[0, 4-rx]  y = 0 to y = 4-rx
    pS = 6*rx*ry; % given
    Qa = int(int(pS, ry, 0, rx), rx, 0, 2);
    Qb = int(int(pS, ry, 0, 4-rx), rx, 2, 4);
    Q = Qa + Qb;
    fprintf('this region/surface has a charge of: %s  C\n', Q);  %  32  C
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    pS = (1e-3)*10*ry*rz*rx^2; % given 
    pS = subs(pS, rz, 4); % since the wedge is at z = 4  ... rz = cz
    
    % rx[0,2] ry[0, (3/2)x ]
    y_ul = (3/2)*rx;
    Q = int(int(pS, ry, 0, y_ul), rx, 0, 2);
    fprintf('charge on triangle is: %.3f  C\n', Q); % 0.288  C
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    pV = (1e-9)*4*cz*cos(cf)*cr^2; % given
    % get charge on the wedge cr[0,2] , cf[0,pi/4] , cz[0,1]
    % use vol:  cr d{cr} d{cf} d{cz}
    intg = cr * pV;
    Q = int(int(int(intg, cr, 0, 2), cf, 0, sym(pi)/4), cz, 0, 1);
    fprintf('charge on wedge volume: %.3f   nC\n', (1e9)*Q);          % 5.657   nC
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    pL = (1e-9)*12*rx^2; % given
    Q = int(pL, rx, 0, 1);
    fprintf('this line has a charge of: %.3f  nC\n', (1e9)*double(Q));  %  4.000  nC
    % r >>> x , so approx
    E = double([ 0, 0, Q / ( 4 * sym(pi) * const.ep0 * 1000^2) ]);
    fprintf('E @ [0, 0, 1000]  about:  [ %.2f , %.2f , %.2f ]  mV/m\n', (1e6).*E);
    % E = [ 0.00 , 0.00 , 36.00 ]  mV/m
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    pL = (1e-6)*12; % given const linear charge density
    % with z = 0 and loop by x^2 + y^2 = 9, --> cr = 3 , cf[0,2pi] , cz = 0
    % only z component effects point charge
    q_pt = (1e-3)*4;
    rq = [0, 0, 4]; % it is 4m above xy plane  so R = 5
    E = double([0, 0, (pL * 3 * 4) / (2 * const.ep0 * 5^3)]); % net E on point charge
    F = double(q_pt .* E);
    fprintf('E at pt:  [ %.3f, %.3f , %.3f]  kV/m\n', (1e-3).*E);
    fprintf('F on pt charge , q E  =  [ %.2f, %.2f, %.2f ]  N\n', F);
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    pS = (1e-9)*12; % given
    % on top of sheet, normal is [ 0 , 0 , 1 ]   reverse for bottom
    norm_top = [0,0,1];
    norm_bot = [0,0,-1];
    E = pS / ( 2 * const.ep0);
    E_top = double(E .* norm_top);
    E_bot = double(E .* norm_bot);
    fprintf('top is E = [ %.3f , %.3f , %.3f ]  V/m\n', E_top);
    fprintf('bot is E = [ %.3f , %.3f , %.3f ]  V/m\n', E_bot);
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    pS = (1e-9)*6; % given
    plane = rx + 2*ry - 5; % f(x,y) provides equation for plane
    poi = [-1, 0, 1];
    gradP = -ee.getGradRec(plane);  % could go +/-  but you want the negative here
    norm_u = gradP ./ norm(gradP);  % turn it into a unit normal vector
    % now it is just applying infinite sheet to normal direction
    E = (pS / (2 * const.ep0)) .* norm_u;
    fprintf('E at poi = [ %.3f , %.3f , %.3f ]  V/m\n', E);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    pS1 = (1e-6)*10; % @ x = 2      E pushes 
    pS2 = (1e-6)*-20; % @ y = -3    E pulls   q < 0
    pS3 = (1e-6)*30; % @ z = 5      E pushes
    
    poi_P = [5, -1, 4];
    r1 = [2, -1, 4]; % closest point on x = 2
    r2 = [5, -3, 4]; % closest point on y = -3
    r3 = [5, -1, 5]; % closest point on z = 5
    r1_P = (poi_P - r1) ./ norm(poi_P - r1); % unit normal x = 2 closest point to poi
    r2_P = (poi_P - r2) ./ norm(poi_P - r2); % unit normal y = -3 closest point to poi
    r3_P = (poi_P - r3) ./ norm(poi_P - r3); % unit normal z = 5 closest point to poi
    E1 = (pS1 / ( 2 * const.ep0)) .* r1_P; % E from x = 2
    E2 = (pS2 / ( 2 * const.ep0)) .* r2_P; % E from y = -3
    E3 = (pS3 / ( 2 * const.ep0)) .* r3_P; % E from z = 5
    E = E1 + E2 + E3;
    fprintf('E at poi P = [ %.1f , %.1f , %.1f ]  kV/m\n', E./(1e3));
    
    poi_R = [0, -2, 1];
    r1 = [2, -2, 1]; % closest point on x = 2
    r2 = [0, -3, 1]; % closest point on y = -3
    r3 = [0, -2, 5]; % closest point on z = 5
    r1_R = (poi_R - r1) ./ norm(poi_R - r1); % unit normal x = 2 closest point to poi
    r2_R = (poi_R - r2) ./ norm(poi_R - r2); % unit normal y = -3 closest point to poi
    r3_R = (poi_R - r3) ./ norm(poi_R - r3); % unit normal z = 5 closest point to poi
    E1 = (pS1 / ( 2 * const.ep0)) .* r1_R; % E from x = 2
    E2 = (pS2 / ( 2 * const.ep0)) .* r2_R; % E from y = -3
    E3 = (pS3 / ( 2 * const.ep0)) .* r3_R; % E from z = 5
    E = E1 + E2 + E3;
    fprintf('E at poi R = [ %.1f , %.1f , %.1f ]  kV/m\n', E./(1e3));
    
    poi_Q = [3, -4, 10];
    r1 = [2, -4, 10]; % closest point on x = 2
    r2 = [3, -3, 10]; % closest point on y = -3
    r3 = [3, -4, 5]; % closest point on z = 5
    r1_Q = (poi_Q - r1) ./ norm(poi_Q - r1); % unit normal x = 2 closest point to poi
    r2_Q = (poi_Q - r2) ./ norm(poi_Q - r2); % unit normal y = -3 closest point to poi
    r3_Q = (poi_Q - r3) ./ norm(poi_Q - r3); % unit normal z = 5 closest point to poi
    E1 = (pS1 / ( 2 * const.ep0)) .* r1_Q; % E from x = 2
    E2 = (pS2 / ( 2 * const.ep0)) .* r2_Q; % E from y = -3
    E3 = (pS3 / ( 2 * const.ep0)) .* r3_Q; % E from z = 5
    E = E1 + E2 + E3;
    fprintf('E at poi Q = [ %.1f , %.1f , %.1f ]  kV/m\n', E./(1e3));
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    pL = (1e-6)*5; % given
    Q = pL * 2 * sym(pi) * 2; % circumfrance is 2 pi r * linear density
    
    % the radius is 2 on yz plane, for poi @ (3, 0, 0)   E pushes [1,0,0]
    push = [3,0,0]; % the ring pushes towards positive x axis  
    R = sqrt(3^2 + 2^2); % from making triangle with ring and points perpendicular distance
    E = (Q/(4 * sym(pi) * const.ep0 * R^3)) .* push;
    D = double(const.ep0 .* E);
    fprintf('D @ (3, 0, 0)  =  [ %.3f , %.3f , %.3f ]   uC/m^2\n', D*(1e6));
    
    % now 2 identical charges Q @ (0, -3, 0) and (0, 3, 0) in addition to ring
    % find Q so D @ P(3,0,0) = [0,0,0]    ...they will have to pull since ring pushes (negative)
    Q = double(-D(1)*4*sym(pi)*18/sqrt(2));
    fprintf('Q must be %.3f   uC\n', Q*(1e6));
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    D = (1e-9) .* [ry^2, 2*rx*ry, -4*rz]; % given electric flux density in nC/m^2
    pV = ee.getDivgRec(D);
    fprintf('the volume charge density is %s  nC/m^3\n', pV*(1e9));
    
    % flux through x = 3 , y[0,6] , z[0,5]
    D = subs(D, rx, 3);
    dS = [1, 0, 0];
    intg = dot(D, dS);
    flux = int(int(intg, ry, 0, 6), rz, 0, 5);
    fprintf('the flux through surface is   %.3f  nC/m^2\n', (1e9)*flux);
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    E = [12*cr*cz*cos(cf), -6*cr*cz*sin(cf), 6*cos(cf)*cr^2];  % given
    D = const.ep0 .* E; % convert to D
    
    % cr[0,2] cf = pi/2 , cz[0,5]   ...a rectangle
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    D = subs(D, cf, sym(pi)/2);
    dS = [0,1,0];
    intg = dot(D, dS);
    flux = double(int(int(intg, cr, 0, 2), cz, 0, 5));
    fprintf('flux across surface:  %.3f   nC\n', flux*(1e9));
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    D = (1e-9).*[ sin(st)*sin(sf), cos(st)*sin(sf), cos(sf)]; % given electric flux density C/m^2
    pV = ee.getDivgSph(D)
    pnt = [2, sym(pi)/6, sym(pi)/3];
    Dpt = subs(pV, [sr, st, sf], pnt);
    fprintf('at point, pV =   %.3f  nC/m^3\n', (1e9).*Dpt);  %0.000  nC/m^3
    
    intg = pV * sin(st) * sr^2;
    flux = int(int(int(intg, sr, 0, 2), st, 0, sym(pi)/6), sf, 0, sym(pi)/3);
    fprintf('flux through surface:  %.5f   nC\n', (1e12)*flux);
    
    %{
    % wants D @ (2, pi/6, pi/3)
    pnt = [2, sym(pi)/6, sym(pi)/3];
    Dpt = subs(D, [sr, st, sf], pnt);
    fprintf('at point, D = [ %.3f , %.3f , %.3f ]  nC/m^2\n', (1e9).*Dpt);
    
    % flux through  sr = 2, st[0,pi/6] , sf[0, pi/3]
    % sA:      <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    >
    dS = [sin(st)*sr^2, 0, 0];
    intg = dot(D, dS);
    intg = subs(intg, sr, 2);
    flux = int(int(intg, st, 0, sym(pi)/6), sf, 0, sym(pi)/3);
    fprintf('flux through surface:  %.3f   nC\n', (1e9)*flux);
    %}
end
    

%------------------------------------------------------------------------------------------ #25
if sel == 25
    D = [8*rx*ry, 4*rx^2, 0];
    pV = simplify(ee.getDivgRec(D));
    fprintf('pV = %s  C/m^3\n', pV);
    
    D = [4*cr*sin(cf), 2*cr*cos(cf), 2*cz^2];
    pV = simplify(ee.getDivgCyn(D));
    fprintf('pV = %s  C/m^3\n', pV);
    
    D = [2*cos(st)/sr^3, sin(st)/sr^3, 0];
    pV = simplify(ee.getDivgSph(D));
    fprintf('pV = %s  C/m^3\n', pV);
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    pV = (1e-3)*12*rx*ry*rz;
    % cube rx[0,2] , ry[0,2] , rz[0,2]
    Qenc = int(int(int(pV, rx, 0, 2), ry, 0, 2), rz, 0, 2);
    fprintf('gauss says charge of cube is: %.3f mC\n', (1e3)*Qenc);
    fprintf(' he also said the enclosed charge = total outward flux\n');
end


%------------------------------------------------------------------------------------------ #27
if sel == 27
    pV1 = 8e-9;  % for sphere r = 1   C/m^2
    pV2 = -6e-3; % for sphere r = 2   C/m^2
    
    %dS [ (sr)^2  sin(st) d{st} d{sf} , 0 , 0 ]
    dS = sin(st) * sr^2;  % should turn into surface area of a sphere
    
    intg = pV1 * dS;
    intg = subs(intg, sr, 1);
    q1 = int(int(intg, st, 0, sym(pi)), sf, 0, 2*sym(pi));
    
    intg = pV2 * dS;
    intg = subs(intg, sr, 2);
    q2 = int(int(intg, st, 0, sym(pi)), sf, 0, 2*sym(pi));
    
    Q = q1 + q2;
    D3 = Q / (4 * sym(pi) * 3^2);
    fprintf('D @ r = 3    [ %.3f , 0 , 0 ]  mC/m^2\n', (1e3)*D3);
    
    %{
    % didn't say C/m^3
    q1 = int(int(int(pV1 * sin(st) * sr^2, sr, 0, 1), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    q2 = int(int(int(pV2 * sin(st) * sr^2, sr, 0, 2), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    Q = q1 + q2;
    D3 = Q / ( 4 * sym(pi) * 3^2);
    fprintf('D @ r = 3    [ %.3f , 0 , 0 ]  mC\n', (1e3)*D3);
    % D @ r = 3    [ -1.778 , 0 , 0 ]  mC
    %}
end


%------------------------------------------------------------------------------------------ #29
if sel == 29
    D = [ 2*rx*ry, rx^2, 0]; % given
    pV = ee.getDivgRec(D); % because pV = divg(D)
    fprintf('volume charge density pV = %s   C/m^3\n', pV);  % pV = 2*ry   C/m^3
    
    % flux on rx[0,1] , ry = 1 , rz[0,1]  just surface
    dS = [ 0 , 1, 0];
    intg = dot(D, dS);
    intg = subs(intg, ry, 1);
    flux = int(int(intg, rx, 0, 1), ry, 0, 1);
    fprintf('flux on surface:  %.3f  C\n', flux);   % 0.333  C
    
    % Q inside:  rx, ry, rz [0,1]
    Q = int(int(int(pV, rx, 0, 1), ry, 0, 1), rz, 0, 1);
    fprintf('cube has an enclosed charge of:  %.3f   C\n', Q);  % 1.000   C
end


%------------------------------------------------------------------------------------------ #30
if sel == 30
    D = (1e-6).*[ 2*cr*(cz+1)*cos(cf) , -cr*(cz+1)*sin(cf) , cos(cf)*cr^2 ]; % given
    pV = ee.getDivgCyn(D);
    fprintf(' pV = %s   uC/m^2\n', pV*(1e6));   % pV = 3*cos(cf)*(cz + 1)   uC/m^2
    
    % confirm Gauss by volume & surface   cr[0,2] , cf[0,pi/2] , cz[0, 4]
    Q = int(int(int(pV * cr, cr, 0, 2), cf, 0, sym(pi)/2), cz, 0, 4);
    fprintf('enclosed charge is:  %.3f   mC\n', Q*(1e6));  % 72.000   mC
    
    %sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >    long way...
    % top , cz = 4     cr d{cr} d{cf}
    dS = [0,0, cr];
    intg = dot(D, dS);
    intg = subs(intg, cz, 4);
    Qtop = int(int(intg, cr, 0, 2), cf, 0, sym(pi)/2);
    % bottom , cz = 0     -cr d{cr} d{cf}
    dS = [0,0, -cr];
    intg = dot(D, dS);
    intg = subs(intg, cz, 0);
    Qbot = int(int(intg, cr, 0, 2), cf, 0, sym(pi)/2);
    % outside ...cylinder part, cr = 2     cr d{cf} d{cz}
    dS = [cr, 0, 0];
    intg = dot(D, dS);
    intg = subs(intg, cr, 2);
    Qout = int(int(intg, cf, 0, sym(pi)/2), cz, 0, 4);
    % side on y=0 plane   cf = pi/2   d{cr} d{cz}
    dS = [0, 1, 0];
    intg = dot(D, dS);
    intg = subs(intg, cf, sym(pi)/2);
    Qy = int(int(intg, cr, 0, 2), cz, 0, 4);
    % side on x=0 plane   cf = 0     -d{cr} d{cz}
    dS = [0,-1,0];
    intg = dot(D, dS);
    intg = subs(intg, cf, 0);
    Qx = int(int(intg, cr, 0, 2), cz, 0, 4);
    Q = Qtop + Qbot + Qout + Qy + Qx;
    fprintf('\nQtop:  %.3f   uC\n', Qtop*(1e6));
    fprintf('Qbot:  %.3f   uC\n', Qbot*(1e6));
    fprintf('Qout:  %.3f   uC\n', Qout*(1e6));
    fprintf('Qy:    %.3f   uC\n', Qy*(1e6));
    fprintf('Qx:    %.3f   uC\n', Qx*(1e6));
    fprintf('enclosed charge is:  %.3f   mC\n', Q*(1e6));  % 72.000   mC
end


%------------------------------------------------------------------------------------------ #33
if sel == 33
    pV = (1e-3)*10/sr^2; % only valid for  1 < r < 4    0 otherwise
    
    % net flux crossing r = 2 is equivelint to pV over 1:2
    % also, D must be [ 10/r , 0 , 0 ] to produce pV
    % sA: (sr)^2  sin(st) d{st} d{sf}   radial   or just enclose it
    D = (1e-3).*[10/sr,0, 0];  % don't do this, you don't know what other values could have been
    %temp = ee.getDivgSph(D); % equals pV
    
    % sr = 2 ,  st[0,pi] , sf[0,2pi]
    intg = pV * sin(st) * sr^2;
    Qenc_v = int(int(int(intg, sr, 1, 2), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    fprintf('by divg/vol @ r = 2, Qenc =  %.3f  mC\n', Qenc_v*(1e3)); % Qenc =  125.664  mC
    %dS = [sin(st)*sr^2, 0, 0];
    %intg = dot(D, dS);
    %intg = subs(intg, sr, 2);
    %Qenc_s = int(int(intg, st, 0, sym(pi)), sf, 0, 2*sym(pi));
    %fprintf('by surf, Qenc = %.3f  mC\n', Qenc_s*(1e3));
    
    % sr = 6
    intg = pV * sin(st) * sr^2;
    Qenc = int(int(int(intg, sr, 1, 4), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    fprintf('by divg/vol @ r = 6, Qenc =  %.3f  mC\n', Qenc*(1e3)); % Qenc =  376.991  mC
    
    fprintf('\n at r = 1 , D = [ 0 , 0 , 0 ]  mC/m^2 ...no pV, no D    D = Q/4 pi r^2\n');
    temp = Qenc / ( 4 * sym(pi) * 5^2); % simple   D = Q / 4 pi r^2
    fprintf(' at r = 5 , D = [ %.3f , 0 , 0 ]  mC/m^2\n', temp*(1e3));
end


%------------------------------------------------------------------------------------------ #35
if sel == 35
    q1 = (1e-9)*2;
    r1 = [1,0,3];
    q2 = (1e-9)*-4;
    r2 = [-2,1,5];
    poi = [1,-2,3];
    mult = 1 / ( 4 * sym(pi) * const.ep0);
    
    r1_p = poi-r1; % direction vector, q1 to poi
    v1 = mult * q1 / norm(r1_p); % potential at poi w/r to q1
    r2_p = poi-r2; % direction vector, q2 to poi
    v2 = mult * q2 / norm(r2_p); % potential at poi w/r to q2
    V = v1 + v2;
    fprintf('potential at point (w/r q1, q2)  is :  %.3f  V\n', V);
end


%------------------------------------------------------------------------------------------ #36
if sel == 36
    q = 8e-9;
    R = sqrt( .03^2 + 2*.02^2);
    V = (4 * q) / ( 4 * sym(pi) * const.ep0 * R);
    fprintf('4 charges make V =  %.2f   volts\n', V); % V =  6985.03   volts
end


%------------------------------------------------------------------------------------------ #37
if sel == 37
    % loop has radius of 4 meters, total charge is Q = 60 uC
    Q = (1e-6)*60;
    rad = 4;
    mult = 1 / ( 4 * sym(pi) * const.ep0 * rad);
    % 2 charges split
    V2 = 2 * (Q/2) * mult;
    fprintf('2 equally spaced charges make V = %.2f kV\n', V2*(1e-3));
    % 3 charges split
    V3 = 3 * (Q/3) * mult;
    fprintf('3 equally spaced charges make V = %.2f kV\n', V3*(1e-3));
    pL = Q / (8*sym(pi)); % because circum = 2 pi r  ...charge per unit length
    dl = mult * pL * 2*sym(pi) * 4 % if it were to be integrated  
end


%------------------------------------------------------------------------------------------ #38
if sel == 38
    q1 = 1e-3;
    r1 = [0,0,4];
    q2 = -2e-3;
    r2 = [-2,5,1];
    q3 = 3e-3;
    r3 = [3,-4,6];
    mult = 1 / ( 4 * sym(pi) * const.ep0);
    
    % find V/potential @ P (-1,1,2)
    P = [-1, 1, 2];
    r1P = P-r1;
    r2P = P-r2;
    r3P = P-r3;
    V1 = mult * q1 * (1/norm(r1P));
    V2 = mult * q2 * (1/norm(r2P));
    V3 = mult * q3 * (1/norm(r3P));
    Vp = V1 + V2 + V3;
    fprintf('potential at P = %.2f  MV\n', Vp*(1e-6));
    
    % to find Vpq  Q (1,2,3)  find V @ Q and subtract
    Q = [1,2,3];
    r1Q = Q-r1;
    r2Q = Q-r2;
    r3Q = Q-r3;
    V1 = mult * q1 * (1/norm(r1Q));
    V2 = mult * q2 * (1/norm(r2Q));
    V3 = mult * q3 * (1/norm(r3Q));
    Vq = V1 + V2 + V3;
    fprintf('potential at Q = %.2f  MV\n', Vq*(1e-6));
    fprintf('Vpq must be VQ - VP = %.2f  MV\n', (Vq-Vp)*(1e-6));
end


%------------------------------------------------------------------------------------------ #39
if sel == 39
    V = sin(cf)*exp(-cz)*cr^2;
    E = -1 * ee.getGradCyn(V);
    poi = [4, sym(pi)/4, -1];
    Ept = subs(E, [cr, cf, cz], poi);
    fprintf('E at point = [ %.2f , %.2f , %.2f ]  V/m\n', Ept);
end


%------------------------------------------------------------------------------------------ #40
if sel == 40
    V = ry*(rz+3)*rx^2;
    E = -1 * ee.getGradRec(V);
    poi = [3, 4, -6];
    Ept = subs(E, [rx, ry, rz], poi);
    fprintf('E at poi : [ %.2f , %.2f , %.2f ]  V/m\n', Ept);
    pV = const.ep0 * ee.getDivgRec(E);
    % want Qenc  x[0,1] , y[0,1] , z[0,1]  ....use divg theorm...no use hitting all surfaces
    Qenc = int(int(int(pV, rx, 0, 1), ry, 0, 1), rz, 0, 1);
    fprintf('the cube has a charge of %.3f  pC\n', Qenc*(1e12));
end


%------------------------------------------------------------------------------------------ #43
if sel == 43
    D = [ 4*rx, -10*ry^2, rz^2 ];
    p = [1, 2, 3];
    divgD = ee.getDivgRec(D);  % maxwell 1,  pV = divg(D)
    pV = subs(divgD, [rx, ry, rz], p);
    fprintf('charge density must be %.2f   C/m^3\n', pV);
end


%------------------------------------------------------------------------------------------ #44
if sel == 44
    q1 = 10e-9;
    r1 = .03;
    area1 = 4*sym(pi)*r1^2;
    D1 = q1/area1;
    
    q2 = -5e-9;
    r2 = .05;
    area2 = 4*sym(pi)*r2^2;
    D2 = (q1+q2)/area2;
    
    fprintf(' r < 3cm, D = [ 0, 0, 0 ]   ...nothing enclosed\n');
    fprintf(' 3cm < r < 5cm  D = [ 0, 0, 10/4 pi r^2 ]  nC/m^3\n');
    fprintf(' r > 5cm  D = [ 0, 0, 5/4 pi r^2 ]  nC/m^3\n');   
end


%------------------------------------------------------------------------------------------ #45
if sel == 45
    syms E0;
    syms a;
    E = [ E0 * cr/a, 0, 0 ]; % for  0 < rho < a
    D = const.ep0 .* E;
    pV = ee.getDivgCyn(D);
end


%------------------------------------------------------------------------------------------ #46
if sel == 46
    E = [ 2*rx*ry*rz, rz*rx^2, ry*rx^2 ];
    V = -1 * potential(E, [rx, ry, rz]);  % big $$$
    start = [ 2, 1, -1 ];
    stop = [ 5, 1, 2 ];
    q = 2e-6;
    Va = subs(V, [rx, ry, rz], start);
    Vb = subs(V, [rx, ry, rz], stop);
    Vab = Vb - Va;
    work = Vab * q;
    fprintf(' work or change in energy is:  %.3f  uJ \n', work*(1e6));
    % or just take it on the path
    path = ee.t01seg(start, stop);
    path_d = diff(path, pt, 1);
    temp = subs(E, [rx, ry, rz], path);
    intg = dot(temp, path_d);
    work = -q*int(intg, pt, 0, 1);
    fprintf(' work or change in energy is:  %.3f  uJ \n', work*(1e6));
    % or be a complete savage
    intg = subs(E(1), [ry, rz], [1,-1]);
    wx = int(intg, rx, 2, 5);
    intg = subs(E(2), [rx, rz], [5,-1]);
    wy = int(intg, ry, 1,1);
    intg = subs(E(3), [rx, ry], [5, 1]);
    wz = int(intg, rz, -1, 2);
    work = -q*(wx+wy+wz);
    fprintf(' work or change in energy is:  %.3f  uJ \n', work*(1e6));
end


%------------------------------------------------------------------------------------------ #47
if sel == 47
    E = [ 12*cr*cz*cos(cf), -6*cr*cz*sin(cf), cos(cf)*6*cr^2 ];
    pV = ee.getDivgCyn(const.ep0 .* E);
    q = 10e-6;
    a = [ 2, sym(pi), -1]; % rec [-2, 0, -1]
    b = [ 2, 0, -1];       % rec [ 2, 0, -1]
    pVa = subs(pV, [cr, cf, cz], a);
    fprintf('charge density @ a = %.3f  nC/m^3\n', pVa*(1e9));
    intg = subs(cr*E(2), [cr, cz], [2,-1]);
    work = -q*int(intg, cf, 0, sym(pi));
    fprintf('work is %.3f  uj\n', work*(1e6));
    
    % should match rectangular from the start
    E = ee.transVecCR(E);
    pV = ee.getDivgRec(const.ep0 .* E);
    a = [-2, 0, -1];
    b = [2, 0, -1];
    pVa = subs(pV, [rx, ry, rz], a);
    fprintf('charge density @ a = %.3f  nC/m^3\n', pVa*(1e9));
    intg = subs(E(1), [ry, rz], [0, -1]);
    work = -q*int(intg, rx, -2, 2);
    fprintf('work is %.3f  uj\n', work*(1e6));   
end

%------------------------------------------------------------------------------------------ #48
if sel == 48
    E = [ 20*sr*sin(st), 10*sr*cos(st), 0];
    q = 10e-9;
    a = [ 5, sym(pi)/6, 0];
    b = [ 5, sym(pi)/2, 0];
    c = [10, sym(pi)/6, 0];
    d = [5, sym(pi)/6, sym(pi)/3];
    e = [10, sym(pi)/2, sym(pi)/3];
    Q = 10e-9;
    % <  d{sr} ,   sr d{st} , sr sin(st) d{sf}  >  a to b is sr = 5, st[30,90], sf = 0
    intg = subs(sr*E(2), [sr, sf], [5, 0]);
    Wab = -Q * int(intg, st, sym(pi)/6, sym(pi)/2);
    fprintf('work a to b =   %.3f  nJ\n', Wab*(1e9));
    % <  d{sr} ,   sr d{st} , sr sin(st) d{sf}  >  a to c is sr[5,10], st=30, sf = 0
    intg = subs(E(1), [st,sf], [sym(pi)/6,0]);
    Wac = -Q * int(intg, sr, 5, 10);
    fprintf('work a to c =   %.3f  nJ\n', Wac * (1e9));
    % <  d{sr} ,   sr d{st} , sr sin(st) d{sf}  >  a to d is sr = 5, st=30, sf[0,60]
    intg = subs(sr*sin(st)*E(3), [sr, st], [5, sym(pi)/6]);
    Wad = -Q * int(intg, sf, 0, sym(pi)/3);
    fprintf('work a to d =   %.3f  nJ\n', Wad * (1e9)); % expected
    % <  d{sr} ,   sr d{st} , sr sin(st) d{sf}  >  a to e is  sr[5,10] , st[30,90] , sf[0,60]
    intg_sr = subs(E(1), [st, sf], [sym(pi)/6, 0]);
    Wsr = -Q * int(intg_sr, sr, 5, 10);
    intg_st = subs(sr*E(2), [sr, sf], [10, 0]);
    Wst = -Q * int(intg_st, st, sym(pi)/6, sym(pi)/2);
    intg_sf = subs(sr*sin(st)*E(3), [sr, st], [10, sym(pi)/2]);
    Wsf = -Q * int(intg_sf, sf, 0, sym(pi)/3);
    work = (Wsr + Wst + Wsf)*(1e9);
    fprintf('work a to e =   %.3f  nJ\n', work);
        %Wae = Wac + Wab + Wad; % no need to do all seperate... NO must account for all transitions
        %fprintf('work a to e =   %.3f  nJ\n', Wae * (1e9));
end


%------------------------------------------------------------------------------------------ #49
if sel == 49
    E = [10/sr^2,0,0];
    a = [1, sym(pi)/4, sym(pi)/2];
    b = [5, sym(pi), 0];
    % <  d{sr}  ,   sr d{st}  , sr sin(st) d{sf}  > 
        %mult = 1 / ( 4 * sym(pi) * const.ep0);   you didn't have Q to begin with...can't do
        %Vab = mult * (( 1/b(1)) - (1/a(1)));
        %fprintf('Vab is %.3f  MV\n', Vab*(1e-6));
   intg_sr = subs(E(1), [st, sf], [sym(pi/4), sym(pi)/2]);
   Vsr = -1 * int(intg_sr, sr, 1, 5)
   % ... E has no st or sf....integration is good
end


%------------------------------------------------------------------------------------------ #50
if sel == 50
    E = [ 20*rx, 40*ry, -10*rz ];
    q = 2e-6;
    a = [2, 0, 0]; % cyn, rec: [2,0,0]
    b = [2, sym(pi)/2, 0]; % cyn, rec: [0, 2, 0]    path shouldn't matter
    % from conv to rec
    intg = subs(E(1), [ry,rz], [0,0]);
    Wx = -q * int(intg, rx, 2, 0);
    intg = subs(E(2), [rx, rz], [0,0]);
    Wy = -q * int(intg, ry, 0, 2);
    work = Wx + Wy;
    fprintf('work is   %.3f  uJ\n', work*(1e6));
    % convert E , keep cyn
    E = ee.transVecRC(E);
    intg = subs(cr*E(2), [cr,cz], [2,0]);
    work = -q * int(intg, cf, 0, sym(pi)/2);
    fprintf('work is   %.3f  uJ\n', work*(1e6));
end


%------------------------------------------------------------------------------------------ #51
if sel == 51
    pS = 40e-9; %  nC/m^2   in x = 0 plane...that means it sits on yz axis
    q = 10e-6;
    a = [3, 4, -1];
    b = [1, 2, 6];
    E = [ pS/(2*const.ep0),0,0];
    work = -q*int(E(1), rx, 3,1);
    fprintf('work is:  %.3f   mJ\n', work*(1e3));
end


%------------------------------------------------------------------------------------------ #52
if sel == 52
    Va = 2*rx^2 + 4*ry^2;
    Ea = (-1)*ee.getGradRec(Va);
    fprintf('Ea = [ %s , %s , %s ]\n', Ea);
    Da = const.ep0 .* Ea;
    pVa = ee.getDivgRec(Da);
    fprintf('pVa = %.3f  nC/m^3\n', pVa*(1e9));
    
    Vb = sin(cf)*10*cr^2 + 6*cr*cz;
    Eb = -1*ee.getGradCyn(Vb);
    fprintf('\nEb = [ %s , %s , %s ]\n', Eb);
    Db =  Eb;
    pVb = ee.getDivgCyn(Db);
    fprintf('pVb = %s *ep0 C/m^3\n', simplify(pVb));
    
    Vc = cos(st)*sin(sf)*5*sr^2;
    Ec = -1 .* ee.getGradSph(Vc);
    fprintf('\nEc = [ %s , %s , %s ] V/m \n', Ec);
    pVc = ee.getDivgSph(Ec);
    fprintf('pVc = ep0 *  %s  C/m^3\n', simplify(pVc));
end

%------------------------------------------------------------------------------------------ #53
if sel == 53
    V = cr*exp(-cz)*sin(cf);
    E = -1 .* ee.getGradCyn(V);
    fprintf('E = [ %s, %s, %s ] V/m\n', E);
    curlE = ee.getCurlCyn(E) % therefore conservative
end


%------------------------------------------------------------------------------------------ #54
if sel == 54
    V = (sin(st)*cos(sf))/sr^3;
    E = -1 .* ee.getGradSph(V);
    D = const.ep0 .* E;
    pnt = [1, sym(pi)/6, sym(pi)/3];
    Dpnt = subs(D, [sr, st, sf], pnt);
    fprintf('D @ point = [ %.3f , %.3f , %.3f ] pC/m^2\n', (1e12).*Dpnt);
end


%------------------------------------------------------------------------------------------ #57
if sel == 57
    D = [ 2*cr*sin(cf), -cos(cf)/(2*cr), 0];
    curlD = ee.getCurlCyn(D);
    fprintf('curl D = [ %s , %s , %s ]  ...not 0,0,0 ...not genuine\n', curlD);  % ???
    %sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} > 
    % cr = 1, cf[0,pi/4] , cz[0,1]
    dS = [ cr, 0, 0];
    intg = dot(D, dS);
    intg = subs(intg, cr, 1);
    flux = int(int(intg, cf, 0, sym(pi)/4), cz, 0, 1);
    fprintf('flux across surface:  %.4f  C\n', flux);
end


%------------------------------------------------------------------------------------------ #60
if sel == 60
    % find p given distance and potential
    p = ( 9 * 4 * sym(pi) * const.ep0 * (1e-9)^2 ) / cos(0);
    % p = Qd,  Q = p/d
    Q = p / (1e-9);
    th = atan( (1e-9)/(1e-9) );
    r = ( 2 * (1e-9)^2 )^(1/2);
    V = double(( p * cos(th) ) / ( 4 * sym(pi) * const.ep0 * r^2));
    fprintf('potential at (1nm,1nm) is %.3f  V\n', V);
end


%------------------------------------------------------------------------------------------ #62
if sel == 62
    p = (1e-6).*[2, 6, -4]; % dipole moment  uC/m
    pos = [2,3,-1];
    poi = [4, 0, 1];
    r = poi - pos;
    V = dot(p,r) / ( 4 * sym(pi) * const.ep0 * norm(r)^2);
    fprintf('V = %.3f  kV\n', (1e-3)*V);   
end


%------------------------------------------------------------------------------------------ #64
if sel == 64
    q1 = 40*(1e-9);
    q2 = -50*(1e-9);
    p1 = [0, 0, 1];
    p2 = [2, 0, 0];
    w1 = 0;
    r12 = p2-p1; % position from q1 to q2
    w2 = q2 * ( q1 / ( 4 * sym(pi) * const.ep0 * norm(r12)));
    fprintf(' total work:  %.3f   uJ\n', (w1+w2)*(1e6));
end


%------------------------------------------------------------------------------------------ #65
if sel == 65
    V = 2*(rx^2) + 6*ry^2;
    E = -1 .* ee.getGradRec(V);
    D = const.ep0 .* E;
    intg = (1/2)*dot(E,D);
    We = int(int(int(intg, rx, -1, 1), ry, -1, 1), rz, -1, 1);
    fprintf('energy stored in volume is   %.3f  nJ\n', (1e9)*We);
end


%------------------------------------------------------------------------------------------ #66
if sel == 66
    E = [ 2*sr*sin(st)*cos(sf), sr*cos(st)*cos(sf), -sr*sin(sf)];
    D = const.ep0 .* E;
    intg = (1/2)*dot(D,E)*sin(st)*sr^2;
    We = int(int(int(intg, sr, 0, 2), st, 0, sym(pi)), sf, 0, sym(pi));
    fprintf('energy in volume is %.3f   nJ\n', (1e9)*We);
end


%------------------------------------------------------------------------------------------ #68
if sel == 68
    E = [ ry^2, 2*rx*ry, -4*rz ];
    D = const.ep0 .* E;
    intg = (1/2)*dot(D, E);
    energy = int(int(int(intg, rx, 0, 2), ry, -1, 1), rz, 0, 4);
    fprintf(' energy in this volume is   %.3f  nJ\n', (1e9)*energy);
end


%------------------------------------------------------------------------------------------ #69
if sel == 69
    V = cr*exp(-cz)*sin(cf);
    E = -1 .* ee.getGradCyn(V);
    D = const.ep0 .* E;
    intg = (1/2) * dot(D, E) * cr;
    enr = int(int(int(intg, cr, 0, 1), cf, 0, sym(pi)*2), cz, 0, 2);
    fprintf(' energy in this volume is   %.3f  pJ\n', (1e12)*enr);
end
    
    
    
    
    
    
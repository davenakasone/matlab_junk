%{
         ch4 notes
      
    #1      ex4.1       coulomb's law
    #2      pp4.1       " "
    #3      ex4.3       use physics to make ODE
    #4      ex4.5       sheet charge
    #5      pp4.5       another square sheet, Q, E, F
    #6      ex4.6       infinite surface and line
    #7      pp4.6       infinite line, considering total radius
    #8      ex4.7       flux density D
    #9      pp4.7       D
    #10     ex4.8       gauss in cyn
    #11     pp4.8       gaus in rec
    #12     pp4.9       good one for spheres inside and outside
    #13     ex4.10      point charges , potential difference
    #14     pp4.10      " "
    #15     ex4.11      potential with point and line charges
    #16     pp4.11      working with potential differences
    #17     ex4.12      given potentials as V scalar, solve for E, D, ect   sph
    #18     pp4.12      "  "  rec
    #19     ex4.13      dipole and potential
    #20     pp4.13      dipole, potential и E
    #21     ex4.14      Energy in a system of 3 point charges
    #22     pp4.14      each charge is positioned one at a time
    #23     pp4.15      energy in a volume
    #24     r2          plane and E
    #25     r5          work
    #26     r10         energy

%}

clc;
close all;
clearvars;


                sel = 26;  % CHANGE CHANGE CHANGE


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
if sel == 1      % find force q1 and q2 put on q3, also the E  at r3
    q1 = 1e-3;
    r1 = [3, 2, -1];
    q2 = -2e-3;
    r2 = [-1, -1, 4];
    q3 = 10e-9;
    r3 = [0, 3, 1];

    F31 = ( (q3*q1) .* (r3-r1) ) ./ ( 4 * pi * const.ep0 * norm(r3-r1)^3 );
    F32 = ( (q3*q2) .* (r3-r2) ) ./ ( 4 * pi * const.ep0 * norm(r3-r2)^3 );
    Fnet = double(F31 + F32);
    fprintf(' net force: [ %d , %d , %d ] mN\n', Fnet*1000);
    fprintf(' E: [ %d , %d , %d ] kV/m\n', Fnet / (q3*1000));
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    q1 = 5e-9;
    r1 = [2, 0, 4];
    q2 = -2e-9;
    r2 = [-3, 0, 5];
    q3 = 1e-9;
    r3 = [1, -3, 7];
    
    F31 = ( (q3*q1) .* (r3-r1) ) ./ ( 4 * pi * const.ep0 * norm(r3-r1)^3 );
    F32 = ( (q3*q2) .* (r3-r2) ) ./ ( 4 * pi * const.ep0 * norm(r3-r2)^3 );
    Fnet = double(F31 + F32);
    fprintf('net force: [ %d , %d , %d ] nN\n', Fnet * (1e9));
    fprintf(' E: [ %d , %d , %d ] V/m\n', Fnet / q3);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % separation with an electric field
    ht = .8; % distance fallen  no initial velocity
    E = 500e3; % 500 kV/m applied
    Q = 9e-6; % charge per meter on all particles  per kg
    syms t;

    % ignored force between particles?
    % differential eqns in terms of time...0 starts imply 0 const
    
    xDisp = @(t) (Q * E / 2)*t.^2;
    yDisp = @(t) -(1/2)*const.g * t.^t;
    
    time = sqrt((.8 * 2)/const.g);
    xd = xDisp(time);
    total = 2 * xd  % meters
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    % rx[0,1] ry[0,1]  rz = 0   
    pS = (1e-9)*rx*ry*( rx^2 + ry^2 + 25 )^(3/2); % given surface density
    
    tempQ = int(pS, rx, 0, 1); % only need in z direction ...it is on xy plane
    Q = int(tempQ, ry, 0, 1); % got the charge
    fprintf(' charge on surface: %.3f  nC\n', Q*(1e9));
    
    % at (0,0,5) the distance h = 5
    % E = Q / R^2 ... still in z direction
    pt = [ 0, 0, 5]; % position vector to point of interest
    rS = [ rx, ry, 0];  % anywhere on sheet  r'
    rd = pt - rS; % establish vector from source (sheet) to point of interest
    intgN = pS * rd;
    intgD = 4 * pi * const.ep0 * ( norm(rd) )^3;
    intg = intgN / intgD; % this is a vector, corresponds to E [ rx, ry, rz ]
    Erx = int(int(intg(1), rx, 0, 1), ry, 0, 1);
    Ery = int(int(intg(2), rx, 0, 1), ry, 0, 1);
    Erz = int(int(intg(3), rx, 0, 1), ry, 0, 1);
    E = [ Erx, Ery, Erz ];
    fprintf(' E = [ %.3f , %.3f , %.3f ]  V/m\n', E);  % E at point (0,0,5)
    
    % the force is just F = qE   q given at -1mC
    fprintf(' force is qE = [ %.5f , %.5f , %.5f ] N\n', ( (-1e-3) ) .* E);  
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    % sheet on rx[-2,2] , ry[-2,2] , z = 0
    pS = 12 * abs(ry) * (1e-3); % given surface charge density
    
    Q = int(int(pS, rx, -2, 2), ry, -2, 2);
    fprintf('charge of surface: %.3f  C\n', Q );
    
    pt = [0, 0, 10];  % position vector to point of interest (point implied)
    rS = [rx, ry, 0]; % position vector to surface 
    rd = pt - rS;     % source to point vector
    intg = ( pS .* rd ) ./ ( 4 * pi * const.ep0 * norm(rd)^3);
    Erx = int(int(intg(1), rx, -2, 2), ry, -2, 2);
    Ery = int(int(intg(2), rx, -2, 2), ry, -2, 2);
    Erz = int(int(intg(3), rx, -2, 2), ry, -2, 2);
    E = [ Erx, Ery, Erz ];
    fprintf('electric field intesity at point is: [ %.3f , %.3f , %.3f ]  MV/m\n', E ./ (1e6));
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % x = 2 carries  10nC/m^2   y = -3 carries 15nC/m^2
    % also a line x = 0, z = 2   10pi nC / m
    % super position each as an infinite plane for net on E @ (1,1,-1)
    
    E = [ 0, 0, 0 ]; 
    temp1 = (10 * (1e-9)) / (2 * const.ep0);
    E(1) = -1 * temp1; % x = 2   plane can only effect -x component
    
    temp2 = (15 * (1e-9)) / (2 * const.ep0);
    E(2) = temp2;  % y = -3   plane can only effect +y component
    
    pnt = [ 1, 1, -1 ];     % point of interest
    source = [ 0, 1, 2 ];   % closest point located on source
    Rv = pnt - source;      % source --> pnt vector
    R = norm(Rv);            % distance between source and point
    RvU = Rv ./ R;
    temp3 = ( 10 * sym(pi) * (1e-9) ) / ( 2 * sym(pi) * const.ep0 * R );
    comps = temp3 .* RvU;
    E = E + comps;
    
    fprintf('E = [ %s , %s , %s ] V/M\n', E); % [ -162 pi , 270 pi , -54 pi ]
end
    
%------------------------------------------------------------------------------------------ #7
if sel == 7    
    % infinite line, x = 0, y = 2     want E @ (1, 1, -1)   pL = 10pi nC / m
    
    pL = 10 * sym(pi) * (1e-9);  % given linear charge density
    pnt = [1, 1, -1];            % point of interest
    source = [0, 2, -1];         % closest source to point of interst      
    rv = pnt - source;           % source to point vector
    R = norm(rv);                % distance between line (source) and point of interest
    rvU = rv ./ R;               % unit vector, source to point
    comp = pL / ( 2 * sym(pi) * const.ep0 * R); % based on infinite line charge
    E = comp .* rvU;
    
    % don't forget to include planes:
    temp1 = (10 * (1e-9)) / (2 * const.ep0);
    E(1) = E(1) + -1 * temp1; % x = 2   plane can only effect -x component
    
    temp2 = (15 * (1e-9)) / (2 * const.ep0);
    E(2) = E(2) + temp2;  % y = -3   plane can only effect +y component
    
    fprintf('E = [ %.3f , %.3f , %.3f ]  V/m\n', E); 
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % want D at (4, 0, 3)       point charge -5pi mC @ (4, 0, 0)
    % and a line charge 3pi mC/m on y axis .... x = 0, z = 0
    
    poi = [ 4, 0, 3 ];  % point of interest...where to calculate D
    
    % start with point charge
    source = [ 4, 0 , 0];  % position of source
    Rv = poi - source;     % getting from source to poi
    R = norm(Rv);          % distance between source and poi
    comp = (-5*sym(pi)*(1e-3)) / ( 4 * sym(pi) * const.ep0 * R^3);
    E1 = comp .* Rv;
    
    % now do line    
    line = [ 0, 0, 0]; % closest point on line to poi
    Rv = poi - line;   % source to poi vector
    R = norm(Rv);      % distance between
    rvU = Rv ./ R;     % unit vector, source to poi
    comp = ( 3*sym(pi)*(1e-3) ) / (2 * sym(pi) * const.ep0 * R);
    E2 = comp .* rvU;
    
    E = E1 + E2;
    D = const.ep0 .* E;
    fprintf('D = [ %.1f , %.1f , %.1f ]    uC / m^2 \n', D .* (1e6) );
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % pt charge 30nC @ origin,  y=3 has 10 nC/m^2   want D @ (0,4,3)
    %burner();
    % handle point charge first
    poi = [ 0, 4, 3];
    source = [0, 0, 0];
    rv = poi - source;
    R = norm(rv);
    rvU = rv ./ R;
    comp = ( 30 * (1e-9) ) / ( 4 * pi * const.ep0 * R^2 );
    E1 = double(comp .* rvU);
    
    % handle infinite plane
    source = [0, 3, 3];  % closest point on source to poi
    rv = poi - source;
    R = norm(rv);
    rvU = rv ./ R;
    comp = ( 10 * (1e-9) ) / ( 2 * const.ep0 );
    E2 = double( comp .* rvU);
    
    E = E1 + E2;
    D = const.ep0 .* E;
    fprintf('flux density D = [ %.3f , %.3f , %.4f ]   nC / m^2\n', D*(1e9) );
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    % given D, want charge density at a point, and charge enclosed cr = 1, cf[0,2pi] cz[-2,2]
    D = [ 0, 0, cz * cr * (cos(cf))^2 ];
    poi = [ 1, sym(pi)/4 , 3];
    fprintf('D at point: [ %s , %s , %s ]    C / m^2\n', subs(D, [cr, cf, cz], poi) );
    
    %sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    dS = [ cr, 0, 0];
    intg = dot(D, dS);
    intg = subs(intg, cr, 1); % looks like 0 , no use continuing  nothing on sides
    
    dS = [0, 0, cr];
    intg = dot(D, dS);
    intg = subs(intg, cz, 2);
    Qtop = int(int(intg, cf, 0, 2*sym(pi)), cr, 0, 1); % top says (2*pi)/3
    
    dS = [0, 0, -cr]; % orientate
    intg = dot(D, dS);
    intg = subs(intg, cz, -2);
    Qbot = int(int(intg, cf, 0, 2*sym(pi)), cr, 0, 1); % top says (2*pi)/3
    
    intg = ee.getDivgCyn(D) * cr;
    Qenc = int(int(int(intg, cr, 0, 1), cf, 0, 2*sym(pi)), cz, -2, 2); % (4*pi)/3 
    
    fprintf('gauss law mathches surface or vol intgl = %s  C \n', Qenc);
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    % find the volume charge density at pt, flux on cube rx[0,1] , ry[0,1] , rz[0,1] 
    % + total charge enclosed  ...it is the flux
    D = [ rz + 2*ry^2 , 4*rx*ry, rx]; % C/m^2 given
    pt = [-1, 0, 3]; % given
    fprintf('D at pt: [ %s , %s, %s ]  C/m^2\n', subs(D, [rx, ry, rz], pt));
    dvgD = ee.getDivgRec(D);
    fprintf('\n...just divg(D)  vol charge density at pt: %s  C/m^3 \n',...
        subs(dvgD, [rx, ry, rz], pt));
    Qenc = int(int(int(dvgD, rx, 0, 1), ry, 0, 1), rz, 0, 1);
    fprintf('\nQ_enc = %s  C\n', Qenc);
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    % given volume charge density = 2r  0<r<10  and 0 otherwise   nC/m^3
    
    % find E at r = 2
    intg = 2*sr*(1e-9) * sr^2 * sin(st); % the volume integral with given density
    Qenc = int(int(int(intg, sr, 0, 2), st, 0, sym(pi)), sf, 0, 2*sym(pi));
    Er = Qenc / ( 4 * sym(pi) * const.ep0 * sr^2 );
    Er = subs(Er, sr, 2); % given radius is 2
    E2m = [ 0, 0, Er]; % only in radial direction
    fprintf('at 2 meters, E = [ %s , %s , %s ] ... [ %s , %s , % .3f ]  V/m\n', E2m, E2m);
    
    % find E at r = 12
    intg = 2*sr*(1e-9) * sr^2 * sin(st); % the volume integral with given density
    Qenc = int(int(int(intg, sr, 0, 10), st, 0, sym(pi)), sf, 0, 2*sym(pi)); % r[10,12] = 0
    Er = Qenc / ( 4 * sym(pi) * const.ep0 * sr^2 );
    Er = subs(Er, sr, 12); % given radius is 2
    E2m = [ 0, 0, Er]; % only in radial direction
    fprintf('at 12 meters, E = [ %s , %s , %s ] ... [ %s , %s , % .3f ]   V/m\n', E2m, E2m);   
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    q1 = -4e-6;       % given charge
    r1 = [ 2, -1, 3]; % given location
    q2 = 5e-6;
    r2 = [ 0, 4, -2 ];
    poi = [ 1, 0, 1]; % find potential here
    
    V1 = (1 / (4 * sym(pi) * const.ep0)) * ( q1 / norm(poi - r1));
    V2 = (1 / (4 * sym(pi) * const.ep0)) * ( q2 / norm(poi - r2));
    V = double(V1 + V2);
    fprintf('assuming 0 potential @ oo ,  V = %.3f   kV\n', V/1000);
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    % V(oo) = 0
    q1 = -4e-6;       % given charge
    r1 = [ 2, -1, 3]; % given location
    q2 = 5e-6;
    r2 = [ 0, 4, -2 ];
    q3 = 3e-6;
    r3 = [0,0,0];
    poi = [ -1, 5, 2]; % find potential here
    
    V1 = (1 / (4 * sym(pi) * const.ep0)) * ( q1 / norm(poi - r1));
    V2 = (1 / (4 * sym(pi) * const.ep0)) * ( q2 / norm(poi - r2));
    V3 = (1 / (4 * sym(pi) * const.ep0)) * ( q3 / norm(poi - r3));
    V = double(V1 + V2 + V3);
    fprintf('assuming 0 potential @ oo ,  V = %.3f   kV\n', V/1000);
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    % point charge of 5 nC @ (-3, 4, 0)        line z=1, y=1  2nC/m
    % potential at any point is V  = Vpt + Vlin
    q_pt = 5e-9;
    r_pt = [-3, 4, 0];
    q_lin = 2e-9;
    orgn = [ 0, 0, 0]; % where V = 0
    
    
    % V = 0 @ (0,0,0)   have to subtract V A(5,0,1) 'poi'
    lin_orgn = norm( orgn - [ 0, 1, 1] ); % line is always (x, 1, 1)     sqrt(2) here
    pt_orgn = norm( orgn - r_pt); % point charge is 5 units from origin
    Vlin_orgn = (-q_lin/( 2 * sym(pi) * const.ep0)) * log(lin_orgn);
    Vpt_orgn = double(q_pt / ( 4 * sym(pi) * const.ep0 * pt_orgn));
    poi = [ 5, 0, 1 ];
    lin_A = norm( poi - [ 5, 1, 1]);
    pt_A = norm( poi - r_pt);
    Vlin_A = (-q_lin/( 2 * sym(pi) * const.ep0)) * log(lin_A);
    Vpt_A = double(q_pt / ( 4 * sym(pi) * const.ep0 * pt_A));
    V = double((Vlin_orgn - Vlin_A) + (Vpt_orgn - Vpt_A));
    Voa = 0 - V   %  8.4766 V     potential difference origin to pt A
    
    % V = 100 at B(1,2,1)   what is V at C(-2, 5, 3)
    B = [1, 2, 1];
    C = [-2, 5, 3];
    pt_B = norm( r_pt - B); % dist point charge to B
    pt_C = norm( r_pt - C); % dist point charge to C
    lin_B = norm( [1,1,1] - B); % dist closest point on line to B
    lin_C = norm( [-2,1,1] - C); % dist closest point on line to C
    Vpt_B = q_pt / ( 4 * sym(pi) * const.ep0 * pt_B);
    Vpt_C = q_pt / ( 4 * sym(pi) * const.ep0 * pt_C);
    Vlin_B = ( -q_lin / ( 2 * sym(pi) * const.ep0 ) ) * log(lin_B);
    Vlin_C = ( -q_lin / ( 2 * sym(pi) * const.ep0 ) ) * log(lin_C);
    V = double(( Vpt_C - Vpt_B ) + ( Vlin_C - Vlin_B ));   % potential B to C
    Vbc = 100 + V % don't forget to add consts   49.8250  V
    
    % if V = -5  @ (0,0,0)  what is Vbc  ?
    % don't even need the referece if common ref assumed...
    V_diff_BC = V  % the absolute difference between these 2 points, no matter what is happening
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    r_pt = [0, 0, 0]; % the point charge is at the origin
    q_pt = 5e-9; % point charge is 5 nC
    r_ref = [0, 6, -8];  % here, V = 2
    
    V =@(dist) double(q_pt / ( 4 * sym(pi) * const.ep0 * dist ));
    
    % want potential at A(-3, 2, 6)
    r_A = [-3, 2, 6];
    pt_ref = norm( r_pt - r_ref); % dist point charge to reference pt
    pt_A = norm( r_pt - r_A); % dist point charge to point A
    V_ref_A_diff = V(pt_A) - V(pt_ref); % the abs potential diff, ref-->A
    V_ref_A = V_ref_A_diff + 2;
    fprintf(' w/r ref--> A, potential at pt A =  %.3f  V\n', V_ref_A);
    
    % want potential at B(1, 5, 7)
    r_B = [1, 5, 7];
    pt_ref = norm( r_pt - r_ref); % dist point charge to reference pt
    pt_B = norm( r_pt - r_B);     % dist point charge to point B
    V_ref_B_diff = V(pt_B) - V(pt_ref); % the abs potential diff, ref-->A
    V_ref_B = V_ref_B_diff + 2;
    fprintf(' w/r ref--> B, potential at pt B =  %.3f  V\n', V_ref_B);
    
    % for Vab just subtract Vb - Va...
    fprintf(' Vab , pt A --> pt B = %.3f  V\n', V_ref_B - V_ref_A);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    V = 10 * sin(st) * cos(sf) / sr^2;
    E = -1 * ee.getGradSph(V); % don't forget to multiply by -1
    
    % a)  find flux density D @ (2, pi/2, 0)
        testingVF = ee.transVecSR(E);
        %fun_graphVF(testingVF, [2,0,0], 0);  % just some graphs to see E
    poi = [2, sym(pi)/2, 0];
    D = const.ep0 .* E;           % D  =  ep0 * E
    D_pt = subs(D, [sr, st, sf], poi);
    fprintf('flux density D at pt = [ %.1f , %s , %s ]  pC/m^2\n', (1e12) .* D_pt);
    
    % for the work A to B, use FTC of line integrals
    ptA = [1, sym(pi)/6, 2*sym(pi)/3];
    ptB = [4, sym(pi)/2, sym(pi)/3];
    q = 10e-6; % charge is 10 uC      Vab = w/q
    funV =@(pnt) subs(V, [sr, st, sf], pnt);
    work = q * ( funV(ptB) - funV(ptA) );  % you already took negative
    fprintf(' the work is: %.3f   uJ\n', work * (1e6)); 
    % the painful way...at least see path doesnt matter
    % sr[1, 4] , st[pi/6, pi/2] , sf[2pi/3, pi/3]
    % length:  <  d{sr} ,   sr d{st}   , sr sin(st) d{sf}  >
    % take it so only one var changes at a time = smart
   
    % sr = 1 , st[pi/6, pi/2] , sf=2 pi/3
    dl_1 = [ 0, sr, 0];
    intg_1 = dot(E, dl_1);
    intg_1 = subs(intg_1, [sr, sf], [1, 2*sym(pi)/3]);
    seg1 = int(intg_1, st, sym(pi)/6, sym(pi)/2 );
    
    % sr = 1 st = pi/2, sf[2 pi/3 , pi/3]
    dl_2 = [0, 0, sr * sin(st)];
    intg_2 = dot(E, dl_2);
    intg_2 = subs(intg_2, [sr, st], [1, sym(pi)/2]);
    seg2 = int(intg_2, sf, 2*sym(pi)/3, sym(pi)/3);
    
    % sr[1,4] st = pi/2, sf = pi/3
    dl_3 = [1, 0, 0];
    intg_3 = dot(E, dl_3);
    intg_3 = subs(intg_3, [st, sf], [sym(pi)/2, sym(pi)/3]);
    seg3 = int(intg_3, sr, 1, 4);
    work = -q * ( seg1 + seg2 + seg3); % the negative camps outside...there are two
    fprintf(' the work is: %.3f   uJ\n', work * (1e6));
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    E = 1000 .* [ry + 3*rx^2, rx, 0];
    V = -potential(E, [rx, ry, rz]); %  V = rx^3 + ry*rx  * 1000
    q = -2e-6;
    ptA = [0, 5, 0];
    ptB = [2, -1, 0];
    
    % the parameter way
    param = ee.t01seg(ptA, ptB); % fundamental method
    param_d = diff(param, pt, 1);
    Et = subs(E, [rx, ry, rz], param);
    intg = dot( Et, param_d );
    work = double(-q * int(intg, pt, 0, 1));
    fprintf('work is : %.1f  mJ\n', work * (1e3));
    
    % the potential way   easy way
    Va = subs(V, [rx, ry, rz], ptA);
    Vb = subs(V, [rx, ry, rz], ptB);
    work = q * (Vb - Va);
    fprintf('work is : %.1f  mJ\n', work * (1e3));
    
    % differential length (0,5,0) --> (2,5,0) --> (2, -1, 0)   long way
    % rx[0,2] , ry = 5 , rz = 0
    dl1 = [1, 0, 0];
    intg1 = dot(E, dl1);
    intg1 = subs(intg1, [ry, rz] , [5, 0]);
    seg1 = int(intg1, rx, 0, 2);
    % rx = 2 , ry[5,-1] rz = 0
    dl2 = [ 0, 1, 0];
    intg2 = dot( E, dl2);
    intg2 = subs(intg2, [rx, rz], [2, 0]);
    seg2 = int(intg2, ry, 5, -1);
    work = double(-q * (seg1 + seg2));
    fprintf('work is : %.1f  mJ\n', work * (1e3));
    
    % y = 5 -3x
    params = [ rx, 5-3*rx, 0];
    params_d = diff(params, rx, 1);
    Ex = subs(E, [rx, ry, rz], params);
    intg = dot(Ex, params_d);
    work = -q * int(intg, rx, 0, 2);
    fprintf('work is : %.1f  mJ\n', work * (1e3));  
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    p1 = (1e-9) .* [0, 0, -5];
    r1 =  [0, 0, -2];
    p2 = (1e-9) .* [0, 0, 9];
    r2 = [0, 0, 3];
    % find potential at (0,0,0) of these 2 dipoles
    poi = [0,0,0];
    V1 =  ( dot(p1 , (poi - r1)) ) ./ ( 4 * sym(pi) * const.ep0 * norm(poi-r1)^3);
    V2 =  ( dot(p2 , (poi - r2)) ) ./ ( 4 * sym(pi) * const.ep0 * norm(poi-r2)^3);
    V = double(V1 + V2);
    fprintf('potential at (0,0,0) from the two diploles: %.2f  V\n', V);
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    p = [ 0, 0, 100e-12];
    r = [ 0, 0, 0];
    
    A = [0, 0, 10];
    Va = ( dot( p, (A-r)) ) / ( 4 * sym(pi) * const.ep0 * norm(A-r)^3 );
    fprintf('the potential at this point is %.1f  mV\n', (1e3) * Va);
    ee.feed(A, 'R');
    temp = ee.pts(3,:);
    th = temp(2);
    mult = norm(p) / ( 4 * sym(pi) * const.ep0 * norm(A-r)^3);
    E = mult .* [2 * cos(st), sin(st), 0];
    E = subs(E, st, th);
    fprintf(' E = [ %.1f , %s , %s]  mV / m\n', E*(1e3));
    
    B = [1, pi/3, pi/2];
    ee.feed(B, 'S');
    B_rec = ee.pts(1,:);
    Vb = dot( p, B_rec - r) / ( 4 * sym(pi) * const.ep0 * norm(B_rec - r)^3);
    fprintf('\nthe potential at this point is %.2f  V\n',  Vb);
    th = B(2);
    mult = norm(p) / ( 4 * sym(pi) * const.ep0 * norm(B_rec - r)^3);
    E = mult .* [2 * cos(st), sin(st), 0];
    E = subs(E, st, th);
    fprintf(' E = [ %.1f , %.4f , %s]  V / m\n', E);
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    q1 = -1e-9;
    r1 = [0, 0, 0];
    q2 = 4e-9;
    r2 = [0, 0, 1];
    q3 = 3e-9;
    r3 = [1, 0, 0];
    
    % method 1
    w1 = 0; % first charge enters system, no work needed
    w2 = q2 * ( q1 / norm( r2 - r1) ); % second charge
    w3 = q3 * ( (q1/norm(r3-r1) )  +   (q2/norm(r3-r2) ) );
    W = double(  ( 1 / (4 * sym(pi) * const.ep0) ) * (w1 + w2 + w3) );
    fprintf('energy in system is %.2f   nJ\n', W*(1e9));
    
    % method 2
    w1 = q1 * ( (q2/norm(r1-r2)) + (q3/norm(r1-r3)) ); % V1 is potential from 2 + 3
    w2 = q2 * ( (q1/norm(r2-r1)) + (q3/norm(r2-r3)) ); % V2 is potential from 1 + 3
    w3 = q3 * ( (q1/norm(r3-r1)) + (q2/norm(r3-r2)) ); % V3 is potential from 1 + 2
    W = double( (1/2) * ( 1 / ( 4 * sym(pi) * const.ep0)) * (w1+w2+w3) );
    fprintf('energy in system is %.2f   nJ\n', W*(1e9));
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    q1 = 1e-9;
    r1 = [ 0, 0, 0];
    q2 = -2e-9;
    r2 = [1, 0, 0];
    q3 = 3e-9;
    r3 = [0, 0, -1];
    q4 = -4e-9;
    r4 = [0, 0, 1];
    
    mult = (1/2) * (1/( 4 * sym(pi) * const.ep0));
    w1 = 0;
    fprintf('after q1 is placed, energy is               %.2f     J  ...   no work done\n', w1);
    tempA = q1 * ( q2/norm(r1-r2)); % q1v1
    tempB = q2 * ( q1/norm(r2-r1)); % q2v2
    w2 = mult * (tempA + tempB); 
    fprintf('after q1 and q2 are placed, energy is       %.2f  nJ\n', w2 *(1e9));
    tempA = q1 * ( (q2/norm(r1-r2)) + (q3/norm(r1-r3)) ); % q1v1
    tempB = q2 * ( (q1/norm(r2-r1)) + (q3/norm(r2-r3)) ); % q2v2
    tempC = q3 * ( (q1/norm(r3-r1)) + (q2/norm(r3-r2)) ); % q3v3
    w3 = mult * (tempA + tempB + tempC);
    fprintf('after q1, q2, and q3 are placed, energy is  %.2f  nJ\n', w3*(1e9));
    tempA = q1 * ( (q2/norm(r1-r2)) + (q3/norm(r1-r3)) + (q4/norm(r1-r4)) ); % q1v1
    tempB = q2 * ( (q1/norm(r2-r1)) + (q3/norm(r2-r3)) + (q4/norm(r2-r4)) ); % q2v2
    tempC = q3 * ( (q1/norm(r3-r1)) + (q2/norm(r3-r2)) + (q4/norm(r3-r4)) ); % q3v3
    tempD = q4 * ( (q1/norm(r4-r1)) + (q2/norm(r4-r2)) + (q3/norm(r4-r3)) ); % q4v4
    w4 = mult * (tempA + tempB + tempC + tempD);
    fprintf('after placing all 4 charges, energy is      %.2f  nJ\n', w4 * (1e9) );
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    V = rx - ry + rx*ry + 2*rz;
    E = -1 * ee.getGradRec(V);
    pnt = [1, 2, 3];
    
    Ept = subs(E, [rx, ry, rz], pnt);
    fprintf('E at point:  [ %s , %s , %s ]  V/m\n', Ept);
    
    % energy in cude x[-1,1] , y[-1,1] , z[-1,1]
    intg = (1/2) * const.ep0 * dot(E,E);
    W = int(int(int(intg, rx, -1, 1), ry, -1, 1), rz, -1, 1);
    fprintf('energy in cube is : %.4f   nJ\n', W*(1e9));
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    source = [0, 0, 10];
    poi = [0, 0, 0];
    normal = poi - source;
    pL = 20e-9;
    E = double((pL / ( 2 * const.ep0)) .* normal)
end


%------------------------------------------------------------------------------------------ #25
if sel == 25
    F = [4, -3, 2];
    disp = [10, 2, -7];
    path = ee.t01seg([0,0,0], disp);
    path_d = diff(path, pt, 1);
    intg = dot(F, path_d);
    work = int(intg, pt, 0, 1);
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    q1 = 1e-6;
    r1 = [-2, 1, 5];
    q2 = 4e-6;
    r2 = [1, 3, -1];
    
    mult = (1/2) * (1/( 4 * sym(pi) * const.ep0 ));
    w1 = 0; % no work
    w1 = q1 * ( q2/norm(r1-r2) );  % or do this
    w2 = q2 * ( q1/norm(r2-r1) );
    W = mult * ( w1 + w2 );
    fprintf('energy stored =  %.2f  mJ\n', W*(1e3));
end
    
    
    
    
    
    


    

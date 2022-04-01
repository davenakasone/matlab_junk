%{
    ee330 ch7 magnetic fields

    #1      ex1     magnetic field from current in a str8 conductor
    #2      pp1     H from straight line not nicley set up to poi
    #3      ex2     H, but this time 2 str8 sources...super position
    #4      pp2     H for semi infinite
    #5      pp3     H at some points around a closed loop
    #6      pp4     H for solenoid
    #7      pp5     H from sheet and K
    #8      pp6     H on a toroid
    #9      ex7.7   A potential on a surface
    #10     pp7.7   A at point and on surf

    #11     hw7, 7.8
    #12     hw7, 7.14
    #13     hw7, 7.36
    #14     hw7, num2

%}
clc;
close all;
clearvars;
format short;


                sel = 13;  % CHANGE CHANGE CHANGE
                
                
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
    I = 10;
    poi = [0,0,5];
    % you want to find H due to side 1...side 1 is a straigh line (000)->(0,0,2)
    % alpha1 = 90°, cos(alpha1)=0  ... -cos(alpha2) = 2/5 ...just use formula
    cr = 5; % implied from poi
    a1 = sym(pi)/2;
    a2 = pi-atan(5/2);
    H_fi = double((I/(4*pi*cr))*(cos(a2)-cos(a1)));
    fprintf('H is in the y_hat direction = %.4f  (A/m)\n', H_fi);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    I = 10;
    % you have the source (1,1,0) to (0,0,0)  ... rho = 5, fi_hat is combonation
    cr = 5;
    a1 = sym(pi)/2;
    a2 = acos(sqrt(2)/sqrt(27));
    acl = [-1,0,0]+[0,-1,0]; % the direction of the linear (source)
    acl = acl ./ norm(acl);
    acr = [0,0,1]; % the implied radial direction is in the z_hat direction
    afi = cross(acl,acr); %  cross linear direction and radial to get fi_hat direction
    H = double(( (I/(4*pi*cr))*(cos(a2)-cos(a1)) ) .* afi)
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    poi = [-3,4,0];
    I = 3;
    
    cr = sqrt(poi(1)^2+poi(2)^2); % consider total distance
    a1 = sym(pi)/2;
    a2 = 0;
    al = [0,0,-1];
    ar = [-3,4,0];
    ar = ar ./ norm(ar);
    af = cross(al,ar);
    af = af./norm(af);
    H2 = double(( (I/(4*pi*cr)) * (cos(a2)-cos(a1)) ) .* af);  % could have just gone cylindrical
    
    cr = 4;
    af = [0,0,1];  %  al = [1,0,0] ...  ar = [0,1,0]
    H1 = double( ((I/(4*pi*cr))*(cos(0)-(3/5))).*af ); 
    
    H = H1 + H2
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    
    poi_a = [2,3,0];
    I = 2;
    al = [0,-1,0];  % straight fillament, orgin to -oo on y axis
    ar = [1,0,0];
    af = cross(al,ar);
    cr = 2;
    H = double( ((I/(4*pi*cr))*((3/sqrt(13))+1)).*af)
    
    poi_b = [3,12,-4];
    ee.feed(poi_b,'R');
        %ee.graph();
    cr = 5;
    al = [0,-1,0];
    ar = [3,0,0]+[0,0,-4];
    ar = ar./norm(ar);
    af = cross(al,ar);
    af = af./norm(af);
    source = [0,0,0];
    R = poi_b-source;
    r = R./norm(R);
    a2 = acos(dot(-1*al,r))
    mult = I/(4*pi*cr)*(cos(a2)+1);
    H = mult.*af
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    % ring of r=.05 z = .01
    % dl = < 0, cr d{cf}, 0 >
    % ...use the ex3 derivation
    I = 50e-3;
    
    poi_1 = [0,0,-.01];
    h = -.02;
    cr = .05;
    H = [0,0, ( (I*cr^2)/(2*(cr^2 + h^2)^(3/2)) ) ]
    
    poi_2 = [0,0,.1];
    h = .09;
    cr = .05;
    H = [0,0, ( (I*cr^2)/(2*(cr^2 + h^2)^(3/2)) ) ]
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    N = 2000;
    l = .75;
    a = .05;
    I = .05;
    n = N/l;
    
    Hz = (n*I/2)*(cos(atan(a/l))-cos(pi/2));   % (0,0,0)   at start
    H = [0,0, Hz]
     
    Hz = (n*I/2)*(cos(pi/2)-cos(pi-atan(a/l)));   % (0,0,.75)   at end
    H = [0,0, Hz]
    
    Hz = (n*I/2)*(cos(atan(a/.25))-cos(pi-atan(a/.5)));   % (0,0,.5)   inside
    H = [0,0, Hz]
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    K = [0,0,.05];  % in y=1 plane
    
    poi1 = [0,0,0];
    n1 = poi1 - [0,1,0];
    H1 = .5 .* cross(K,n1)
    
    poi2 = [1,5,-3];
    n2 = poi2 - [1,1,-3];
    n2 = n2 ./ norm(n2);
    H2 = .5 .* cross(K,n2)
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    cent = [0,0,0];
    N = 1000;
    cr0 = .1;
    a = .01;
    I = .1;
    
    
    poi1 = [.03,-.04,0]; % it is not in the coils
    rad = norm(poi1-cent); % this is outside the range [.09,.11]
    H1 = 0;
    
    poi2 = [.06,.09,0];
    rad = norm(poi2-cent); % this is inside: 0.1082
    H2 = (N*I)/(2*pi*rad)
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    %<  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    dS = [0,1,0];        % cr[1,2]   cf = pi/2  cz[0,5]
    A = [0,0,(-cr^2)/4];
    B = ee.getCurlCyn(A);
    intg = dot(B,dS);
    
    flux = int(int(intg, cr, 1,2),cz,0,5);
    fprintf('flux is %.3f  Wb\n', flux);  % or be a savage and do 4 line integrals
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    A = [ry*rx^2, rx*ry^2, -4*rx*ry*rz];
    B = ee.getCurlRec(A);
    poi = [-1,2,5];
    Bpt = subs(B,[rx,ry,rz],poi)
    
    % rx[0,1]  ry[-1,4]  z = 1     ...
    dS = [0,0,1];
    intg = dot(B,dS);
    intg = subs(intg, rz, 1);
    flux = int(int(intg, rx, 0, 1), ry, -1,4)
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    I = 10;
    
    Rmag2 = sqrt((-rx)^2 + (rx-2)^2 + 5^2);
    R2 = [-rx, rx-2, 5];
    dl2 = [1,-1,0];
    cp2 = ee.curlRec2(dl2, R2);
    
    intgx2 = (cp2(1)*I)/(4*sym(pi)*Rmag2^3);
    H2x = double(int( intgx2, rx,2,1));
    H2y = H2x;
    intgz2 = (cp2(3)*I)/(4*sym(pi)*Rmag2^3);
    H2z = double(int(intgz2, rx,2,1));
    H2 = [ H2x, H2y, H2z ];
    
    Rmag1 = sqrt((-rx)^2 + 5^2);
    R1 = [-rx, 0, 5];
    dl1 = [1,0,0];
    cp1 = ee.curlRec2(dl1, R1);
    
    intgy1 = (cp1(2)*I)/(4*sym(pi)*Rmag1^3);
    H1x = 0;
    H1z = 0;
    H1y = double(int(intgy1, rx,0,2)); 
    H1 = [ H1x, H1y, H1z ];
    
    Rmag3 = sqrt(2*(-rx)^2 + 5^2);
    R3 = [-rx, -rx, 5];
    dl3 = [1,1,0];
    cp3 = ee.curlRec2(dl3, R3);
    H3z = 0;
    intgx3 = (cp3(1)*I)/(4*sym(pi)*Rmag3^3);
    H3x = double(int(intgx3, rx,1,0));
    intgy3 = (cp3(2)*I)/(4*sym(pi)*Rmag3^3);
    H3y = double(int(intgy3, rx,1,0));
    H3 = [ H3x, H3y, H3z ];
    
    H = H1 + H2 + H3;
    fprintf('H = [ %.3f , %.3f , %.3f ]  mT\n', H(1)*(1e3), H(2)*(1e3), H(3)*(1e3));
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    I = 10;
    r = .04;
    Rp = [0,0,0];
    
    Rs1t = [ rx, sqrt(r^2 -rx^2), 0 ];
    R1t = Rp - Rs1t;
    dl1t = diff(Rs1t,rx,1);
    cp1t = ee.crossRec(dl1t, R1t);
    intgz1t = (cp1t(3)*I)/(4*sym(pi)*r^3);
    H1t = [0,0,double(int(intgz1t, rx,0,-.04))];
    
    Rs1b = [ rx, -sqrt(r^2 -rx^2), 0 ];
    R1b = Rp - Rs1b;
    dl1b = diff(Rs1b,rx,1);
    cp1b = simplify(ee.crossRec(dl1b, R1b));
    intgz1b = (cp1b(3)*I)/(4*sym(pi)*r^3);
    temp = int(intgz1b,rx);
    H1b = [0,0,double(int(intgz1b, rx,-.04,0))];
    H1 = H1t+H1b;
    
    Rs2 = [ rx, -.04, 0 ];
    R2 = Rp - Rs2;
    R2mag = sqrt( R2(1)^2 + R2(2)^2 + R2(3)^2 );
    dl2 = diff(Rs2,rx,1);
    cp2 = simplify(ee.crossRec(dl2, R2));
    intgz2 = (cp2(3)*I)/(4*sym(pi)*R2mag^3);
    %pretty(intgz2);
    temp = int(intgz2,rx);
    %pretty(temp)
    H2 = [0,0,double(int(intgz2, rx,0,1))];
    
    Rs3 = [ 1, ry, 0 ];
    R3 = Rp - Rs3;
    R3mag = sqrt( R3(1)^2 + R3(2)^2 + R3(3)^2 );
    dl3 = diff(Rs3,ry,1);
    cp3 = simplify(ee.crossRec(dl3, R3));
    intgz3 = (cp3(3)*I)/(4*sym(pi)*R3mag^3);
    %pretty(intgz3);
    temp = int(intgz3,ry);
    %pretty(temp)
    H3 = [0,0,double(int(intgz3, ry,-.04,.04))];
    
    Rs4 = [ rx, .04, 0 ];
    R4 = Rp - Rs4;
    R4mag = sqrt( R4(1)^2 + R4(2)^2 + R4(3)^2 );
    dl4 = diff(Rs4,rx,1);
    cp4 = simplify(ee.crossRec(dl4, R4));
    intgz4 = (cp4(3)*I)/(4*sym(pi)*R4mag^3);
    %pretty(intgz4);
    temp = int(intgz4,rx);
    %pretty(temp);
    H4 = [0,0,double(int(intgz4, rx,1,0))];
    
    H = H1 + H2 + H3 + H4  
end

%------------------------------------------------------------------------------------------ #13
if sel == 13
    cf_ll = 0;
    cf_ul = (5*sym(pi))/18;
    
    %cf_ll = (90-25)*sym(pi)/180;  % does not work...wtf?
    %cf_ul = (90+25)*sym(pi)/180;
    
    cz_ll = 0;
    cz_ul = .2;
    H = [ (10^6)*sin(2*cf)/cr,0,0];
    B = const.mu0 .* H;
    dS = [cr, 1, cr];
    intg = dot(B,dS);
    fz = int(intg, cz, 0, .2);
    flux = double(int(fz, cf, cf_ll, cf_ul))
    
    mult = (-1*sym(pi))/25;
    ll = cos(2*cf_ll);
    ul = cos(2*cf_ul);
    intgl = ul-ll;
    ingll = double(intgl*mult)
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    syms a; assume(a, 'real');
    syms b; assume(b, 'real');
    syms zp; assume(zp, 'real');
    syms I; assume(I, 'real');
    
    Rs = [b+a*cos(cf), cf, 0];
    Rp = [0,0,zp];
    R = Rp-Rs;
    Rmag = sqrt( R(1)^2 + R(3)^2 );
    dl = diff(Rs,cf,1);
    dl(2)=dl(2)*(b+a*cos(cf));
    cp = cross(dl,R);
    
    test1c = [1,0,0];
    test1r = ee.simpleCRsub(test1c)
    test2c = [1,sym(pi)/2,0];
    test2r = ee.simpleCRsub(test2c)
    cp_test = cross(test1r,test2r)
    p_test = ee.crossCyn(test1c,test2c)
end
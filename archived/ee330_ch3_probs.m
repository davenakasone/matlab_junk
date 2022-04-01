%{
    chapter 3 problems, supposed to be the key to the class
        be comfortable handling line integrals by total differential or parameterization
        line integrals on scalar (make curtain) are not too common...vector is the focus

        #1      p3.4r
        #2      p3.7r
        #3      p3.1    lengths
        #4      p3.2    surface areas
        #5      p3.3    volumes
        #6      p3.4    length of path
        #7      p3.5    sa on sph
        #8      p3.6    volume
        #9      p3.7    line integral on parabola
        #10     p3.8    line integral between points   parameterize
        #11     p3.9    line integral, work
        #12     p3.10   line integral on segments
        #13     p3.11   line integral                  displays params vs total diff well
        #14     p3.12   ???
        #15     p3.13   flux on surface rec
        #16     p3.14   flux on volume rec
        #17     p3.15   flux on volume sph
        #18     p3.16   general scalar volumes on rec и cyn
        #19     p3.17   scalar gradients
        #20     p3.18   gradient of scalar eval at pt
        #21     p3.19   a proof
        #22     p3.20   using gradient to find direction of maximum change
        #23     p3.21   use gradient on family of planes....it is normal vector
        #24     p3.22   implication of gradient
        #25     p3.23   directional derivative using gradient as normal
        #26     p3.24   angle between 2 planes...use normals
        #27     p3.25   simple scalar --> gradient @ point
        #28     p3.26   a proof  grad(UV) = U grad(V) + V grad(U)
        #29     p3.27   simple divg
        #30     p3.28   divg @ pt
        #31     p3.29   implication of heat flow vector to scalar temperature
        #32     p3.30   proof of a divg property
        #33     p3.31   nambla operations
        #34     p3.32   flux...is not circulation....don't be trying to stokes for circulation 
        #35     p3.33   divg compare surface and just using theorem  CHECK
        #36     p3.34   more flux on volume
        #37     p3.35   "          "
        #38     p3.36   "          "
        #39     p3.37   required a vector trans  or just be savage and take it over a circle
        #40     p3.38   divg verification...
        #41     p3.39    a hollow object...use superposition on divg, or just do surfaces
        #42     p3.40   curl...
        #43     p3.41   curl
        #44     p3.42   curl with order of operations
        #45     p3.43   curl of gravity vf is 0
        #46     p3.44a  validating stokes theorem vs circulation  rec
        #47     p3.44b  validating stokes theorem vs circulation   cyn
        #48     p3.45   stokes instead of a bunch of line integrals
        #49     p3.46   verify stokes
        #50     p3.47   curl tricks
        #51     p3.48   algebra tricks
        #52     p3.49   a culmination problem
        #53     p3.51   operations
        #54     p3.52   dvg/crl
        #55     p3.53   proofs
        #56     p3.54   proofs
        #57     p3.55   laplacians
        #58     p3.56   laplacian at point
        #59     p3.57   position vector r
        #60     p3.58   laplacians
        #61     p3.59   laps
        #62     p3.60   properties
        #63     p3.61   proof lap = divg(grad)
        #64     p3.62   vector laplcian cyn
        #65     p3.63   property
        #66     p3.64   classify vector feilds
        #67     p3.65   analyzing vector fields
        #68     p3.66   electric field due to line charge solendial and consv


%}
clc;
close all;
clearvars;


                       sel = 68;  % CHANGE CHANGE CHANGE


ee = cls_EE330_helper();
global rx; global ry; global rz; % rectangular params  rx, ry, rz        
global cr; global cf; global cz; % cylindrical params  cr, cf, cz
global sr; global st; global sf; % spherical params    sr, st, sf
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); syms cf; assume(cf,'real'); syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); syms st; assume(st,'real'); syms sf; assume(sf, 'real');
global p_t; syms p_t; assume(p_t, 'real'); % paramater t for a space curve


%------------------------------------------------------------------------------------------ #1
if sel == 1
    vec = [rx, ry, rz];  % given the standard positioning vector
    %fun_graphVF(vec, [0,0,1], 1);
    mag = sqrt(sum(vec.^2));
    gradMag = ee.getGradRec(mag);
    fprintf('the gradient of the magnitude is:  \n[ %s ,\n  %s ,\n  %s ]\n', simplify(gradMag));
    
    divg = ee.getDivgRec(vec);
    fprintf('\nthe divergence = %.2f\n', divg);
    
    dots = dot(vec, vec);
    lapDots = ee.getLaplacianRecS(dots);
    fprintf('\nthe Laplacian of the vector dotted with itself is: %s\n', lapDots);
    
    crl = ee.getCurlRec(vec);
    fprintf('\nthe curl is  [ %.1f , %.1f, %.1f ]\n', crl);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    A = [3*(rx^2)*ry*rz, (rx^3)*rz, (rx^3)*ry-2*rz]; % given
    divgA = ee.getDivgRec(A);
    fprintf('divg = %s         divg not 0, not solenoidal\n', divgA);
    curlA = ee.getCurlRec(A);
    fprintf('curl = [ %.0f, %.0f, %.0f ]   ...  no curl = conservative\n', curlA);
    vecLap = ee.getLaplacianRecV(A);
    fprintf('vector Laplacian: [ %s , %s , %s ]\n', vecLap);
    %fun_graphVF(A, [0,0,1], 1);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    %ee.help_diff();
    
    % cyn length diffs:  <  d{cr}           ,  cr d{cf}     ,  d{cz}          >
    % cr = 3, cf[pi/4, pi/2] , cz = const ... only change on cf
    lenA = int(3, cf, sym(pi)/4, sym(pi)/2);
    fprintf(' the length is: %s\n', lenA);
    
    % sph length diffs:   < d{sr} ,   sr d{st}  , sr sin(st) d{sf}  >
    % sr = 1 , st = pi/6, sf[0, pi/3]
    lenB = int( 1 * sin(pi/6), sf, 0, sym(pi)/3);
    fprintf(' the length is: %s\n', lenB);
    
    % sr = 4, st[pi/6, pi/2], sf = const
    lenC = int(4, st, sym(pi)/6, sym(pi)/2);
    fprintf(' the length is: %s\n', lenC);
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    %ee.help_diff();
    % cyn SA diffs: <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    % cr = 2, cf[pi/3, pi/2], cz[0, 5]  orientation is good 
    temp = int( 2 , cf, sym(pi)/3, sym(pi)/2);
    saA = int(temp, cz, 0, 5);
    fprintf('saA = %s\n', saA);
    
    % cr[1, 3] , cf[0, pi/4] , cz = 1     orientation is bad
    temp = int(cr, cr, 1, 3);
    saB = int(temp, cf, 0, sym(pi)/4);
    fprintf('saB = %s\n', saB);
    
    % sph SA diffs: < (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}  >
    % sr = 10, st[ pi/4, 2pi/3 ]  , sf[ 0 , 2pi ]   orientation good
    temp = int( (10^2)*sin(st), st, sym(pi)/4, 2*sym(pi)/3);
    saC= int(temp, sf, 0, 2*sym(pi));
    fprintf('saC = %s   .... %.4f\n', saC, double(saC));
    
    % sr[0, 4]  ,  st[pi/3, pi/2] , sf = const
    temp = int(sr, sr, 0, 4);
    saD = int(temp, st, sym(pi)/3, sym(pi)/2);
    fprintf('saD = %s     ... %.3f\n', saD, double(saD));
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    %ee.help_diff();
    % rec, vol:  d{rx} d{ry} d{rz}  SCALAR
    % rx[0,1] , ry[1,2], rz[-3,3]
    volA = int(int(int(1, rx, 0, 1), ry, 1, 2), -3, 3);
    fprintf('volA = %.3f\n', volA);
    
    % cyn, vol:  cr d{cr} d{cf} d{cz}  SCALAR
    % cr[2, 5], cf[pi/3, pi], cz[-1,4]
    volB = int(int(int(cr, cr, 2, 5), cf, sym(pi)/3, sym(pi)), cz, -1, 4);
    fprintf('volB = %s     ... %.4f\n', volB, double(volB));
    
    % sph, vol:  (sr)^2 sin(st) d{sr} d{st} d{sf}  SCALAR
    % sr[1, 3], st[pi/2, 2pi/3] , sf[pi/6, pi/2]
    volC = int(int(int( sin(st) * sr^2, sr, 1, 3), st, sym(pi)/2, 2*sym(pi)/3),...
        sf, sym(pi)/6, sym(pi)/2);
    fprintf('volC = %s   ... %.3f\n', volC, double(volC));
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    p1 = [4, 0, 0];
    p2 = [4, pi/6, 0];
    %ee.help_diff();
    % cyn length: <  d{cr}           ,  cr d{cf}     ,  d{cz}          >
    path = int(4, cf, 0, sym(pi)/6); % .... length of  a sector
    fprintf(' point to point is: %s  ...  %.3f\n', path, double(path)); 
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    % sph SA diffs: < (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}  >
    % sr = 5, st[0, pi/4] , sf[0, pi/2]
    temp = int(int( (5^2)*sin(st), st, 0, sym(pi)/4), sf, 0, sym(pi)/2);
    fprintf('sa = %s     .... %.4f\n', temp, double(temp));
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % cyn, vol:  cr d{cr} d{cf} d{cz}  SCALAR
    % cr[2,5], cf[0, pi/3], cz[0,10]
    vol = int(int(int(cr, cr, 2, 5), cf, 0, sym(pi)/3), cz, 0, 10);
    fprintf('vol = %s    ...%.3f\n', vol, double(vol));
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    Hv = [rx*ry^2, ry*rx^2, 0];   % solve on curve rx = ry^2
    %fun_graphVF(Hv, [0,0,1], 1);
    % you should be able to see the gradient --> scalar is f(x,y) = (1/2) rx^2 ry^2
    start = [1, 1, 0];
    stop = [16, 4, 0];
    poten = (1/2)*(rx^2)*(ry^2);
    lower = subs(poten, [rx, ry, rz], start);
    upper = subs(poten, [rx, ry, rz], stop);
    difernce = upper - lower;
    fprintf('the distance = %.3f\n', difernce); % by FTC of line integrals
    %or
    % vector field is transfered with parameters given by path
    vec = subs(Hv, rx, ry^2);
    dl = [2*ry, 1, 0];
    integrand = dot(vec, dl);
    length = int(integrand, 1, 4);
    pot = potential(Hv, [rx, ry, rz]); % big $$$
    fprintf('dist =         %.3f\n', length);
    display(pot);
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    fld = [2*(rx^2)-4*rx*ry, 3*rx*ry-2*ry*rx^2, 0];
    start = [1, -1, 2]; % given
    stop = [3, 1, 2]; % given
    path = stop - start; % implied
        %fun_graphVF(fld, [start;stop], 1);
    fprintf('the curl is not [0, 0, 0] no ftc, not conservative = [%s, %s, %s ]\n', ee.getCurlRec(fld)); 
    % need to make a vector line with parameter t   assume [0, 1]
    rt = (1-p_t).*start + p_t.*stop;
    fprintf('the space curve is represented by [ %s, %s, %s ]  for [rx, ry, rz]\n', rt);
        %fun_graph_spaceC(rt,[0,1]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    fprintf('the differential length is [ %s, %s, %s ] --> [ d{rx}, d{ry}, d{rz} ]\n', dl );
    fld_t = subs(fld, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    fprintf('the vector field converted to param_t : \n[ %s ,\n  %s ,\n  %s ]\n', fld_t);
    integrand = dot(fld_t, dl);
    fprintf('the integrand formed = %s\n', simplify(integrand));
    leng = int(integrand, p_t, 0, 1);
    fprintf('the line integral is    %.3f        only had to do t[0,1]\n', leng);   
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    Fv = [2*rx*ry, (rx^2)-(rz^2), -3*rx*(rz^2)];
    divgF = ee.getDivgRec(Fv);
    %curlF = ee.getDivgRec(Fv);   not 0
    %pot = potential(Fv, [rx, ry, rz]);  failed
    ptA = [0, 0, 0];
    way1 = [0, 1, 0];
    way2 = [2, 1, 0];
    ptB = [2, 1, 3];
        %fun_graphVF(Fv, [ptA;way1;way2;ptB], 1);
        
    % ptA to way1
    rt = (1-p_t).*ptA + p_t.*way1;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    fld_t = subs(Fv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(fld_t, dl);
    leng1 = int(integrand, p_t, 0, 1);
    fprintf('first segment: %.3f\n', leng1);
    
    % way1 to way2
    rt = (1-p_t).*way1 + p_t.*way2;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    fld_t = subs(Fv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(fld_t, dl);
    leng2 = int(integrand, p_t, 0, 1);
    fprintf('second segment: %.3f\n', leng2);
    
    % way2 to ptB
    rt = (1-p_t).*way2 + p_t.*ptB;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    fld_t = subs(Fv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(fld_t, dl);
    leng3 = int(integrand, p_t, 0, 1);
    fprintf('third segment: %.3f\n', leng3);
    fprintf('total is: %.3f\n', leng1 + leng2 + leng3);
    
    % another option, by differentials
    integrand = [ diff(Fv(1), rx, 1), diff(Fv(2), ry, 1), diff(Fv(3), rz, 1) ];
    % ptA to way1    rx[0, 2]    ry[0, 1]   rz[0,3]
    tempX = int(int(integrand(1), ry, 0, 1), rz, 0, 3);
    tempY = int(int(integrand(2), rx, 0, 2), rz, 0, 3);
    tempZ = int(int(integrand(3), rx, 0, 2), ry, 0, 1);
    leng = tempX + tempY + tempZ;
    fprintf('alternatively, = %s\n', leng);
    
    % striaght
    rt = (1-p_t).*ptA + p_t.*ptB;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    fld_t = subs(Fv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(fld_t, dl);
    leng = int(integrand, p_t, 0, 1);
    fprintf('\nstraight path: %.3f  ...work done\n', leng);
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    Fv = [cr^2, cz, cos(cf)];
    ptP = [2, 0, 0];
    ptQ = [2, pi/4, 3];
    % path cr = 2, cf[0, pi/4] , z = 0      --> cr = 2, cf = pi/4, cz[0, 3]
    % cyn length diffs:  <  d{cr} ,  cr d{cf}  ,  d{cz}   >
    
    % frist part cr = 2, cf[0, pi/4] , z = 0        
    dl = [0, cr, 0]; % only movt on acf
    temp = subs(Fv, [cr, cz], [2, 0]);
    integrand = dot(temp, dl);
    part1 = int(integrand, cf, 0, sym(pi)/4);
    fprintf('first part = %s      ...  %.3f\n', part1, double(part1));
    
    % second part cr = 2, cf = pi/4, cz[0, 3]       
    dl = [0, 0, 1]; % only movt on acz
    temp = subs(Fv, [cr, cf], [2, pi/4]);
    integrand = dot(temp, dl);
    part2 = int(integrand, cz, 0, 3);
    fprintf('second part = %s      ...  %.3f\n', part2, double(part2)); 
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    Hv = [rx-ry, rz*ry+rx^2, 5*ry*rz];
    ptA = [1, 0, 0];
    ptB = [0, 0, 0];
    ptC = [0, 0, 1];
    ptD = [0, 2, 0];
        %fun_graphVF(Hv, [ptA;ptB;ptC;ptD], 1);   ...nothing to  z component
    pot = potential(Hv, [rx, ry, rz]); % NaN
    curlH = ee.getCurlRec(Hv); % not 0
    
    % A to B
    rt = (1-p_t).*ptA + p_t.*ptB;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    Ht = subs(Hv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(Ht, dl);
    len1 = int(integrand, p_t, 0, 1);
    fprintf('first seg: %.2f\n', len1);
    
    % B to C
    rt = (1-p_t).*ptB + p_t.*ptC;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    Ht = subs(Hv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(Ht, dl);
    len2 = int(integrand, p_t, 0, 1);
    fprintf('second seg: %.2f\n', len2);
    
    % C to D
    rt = (1-p_t).*ptC + p_t.*ptD;
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)];
    Ht = subs(Hv, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    integrand = dot(Ht, dl);
    len3 = int(integrand, p_t, 0, 1);
    fprintf('second seg: %.2f\n', len3);
    fprintf('\t\t\ttotal : %.2f\n', len1 + len2 + len3);
    
    fprintf('\nalternativley, instead of parameterizing, use total differential\n');
    % or rx[1,0] ry=0=rz     A to B
    temp = Hv(1); % only concerned about x comp in this direction
    integrand = subs(temp, [ry, rz], [0, 0]); % should just be left with rx
    comp1 = int(integrand, rx, 1, 0);   % the limits of integration handle sign change
    fprintf('the rx component: %.2f\n', comp1);
    
    % rx = 0 , ry = 0 , rz[0, 1]   B to C
    temp = Hv(2); % only concerned about y comp in this direction
    integrand = subs(temp, [rx, ry], [0, 0]); % just left with 0....integrating 0 = 0
    comp2 = int(integrand, ry, 0, 1); 
    fprintf('the ry compoenent: %.2f\n', comp2);
    
    % rx = 0, ry[0, 2] , rz[1, 0]  C to D 
    % but notice that    rz = -(1/2) ry + 1    , ry = ( rz - 1) * (-2) = 2 - 2rz    !!!
    tempY = Hv(2);
    tempZ = Hv(3);
    integrandY = subs(tempY, rx, 0); % should be rz*ry
    integrandY = subs(integrandY, rz, -(1/2)*ry + 1); % should be  -ry*(ry/2 - 1)
    integrandZ = subs(tempZ, rx, 0); % should be 5*ry*rz
    integrandZ = subs(integrandZ, ry, 2-2*rz); % should be   -5*rz*(2*rz - 2)
    comp3y = int(integrandY,ry, 0, 2);
    comp3z = int(integrandZ, rz, 1, 0);
    fprintf('composite component for last segment: %.2f\n', comp3y + comp3z);
    fprintf('\t\t\ttotal : %.2f\n', comp1 + comp2 + comp3y + comp3z);
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    Bv = [rx*ry, -ry*rz, rx*rz]; % circullation x = 1...or y=1??? 
    ptA = [0, 0, 0]; % -->
    ptB = [1, 0, 0]; % -->
    ptC = [1, 0, 1]; % --> ptA
    % from stokes, cirrculation --> surf(curl(Bv))
    integrand = ee.getCurlRec(Bv); % but the area is 1/2 base * height = 1/2
    fprintf('circulation must be (1/2) * [ %s, %s, %s ]\n', integrand);
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    Av = [ry, rz, rx]; % want flux on rx[0,1] , ry=1 , rz[0,2]
    % length:  <  d{rx}        ,  d{ry}        ,  d{rz}       >  VECTOR
    % sA:      <  d{ry} d{rz}  ,  d{rx} d{rz}  ,  d{rx} d{ry} >  VECTOR
    % could use stokes or take circulation on each side
    ptA = [1, 1, 0];
    ptB = [0, 1, 0];
    ptC = [0, 1, 2];
    ptD = [1, 1, 2]; %--> ptA
        %fun_graph_path3( [ptA; ptB; ptC; ptD; ptA], Av);
    curlA = ee.getCurlRec(Av);
    % the area is const on ry direction, so use [0, 1, 0] --> d{rx} d{rz}
    dS = [0, 1, 0]; 
    integrand = dot(curlA, dS); % no need to sub y = 1...curl is all const = [-1, -1, -1]
    flux = int(int(integrand, rx, 0, 1), rz, 0, 2);
    fprintf('stokes, flux on surface is: %.1f\n', flux);
    
    % confirm by param line int sum   p_t[0,1]  
    % A to B
    rt = (1-p_t).*ptA + p_t .* ptB;
    fld_t = subs(Av, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)]; 
    integrand = dot(fld_t, dl);
    segAB = int(integrand, p_t, 0, 1);
    fprintf('\nAB segment: %.1f\n', segAB);
    % B to C
    rt = (1-p_t).*ptB + p_t .* ptC;
    fld_t = subs(Av, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)]; 
    integrand = dot(fld_t, dl);
    segBC = int(integrand, p_t, 0, 1);
    fprintf('BC segment: %.1f\n', segBC);
    % C to D
    rt = (1-p_t).*ptC + p_t .* ptD;
    fld_t = subs(Av, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)]; 
    integrand = dot(fld_t, dl);
    segCD = int(integrand, p_t, 0, 1);
    fprintf('CD segment: %.1f\n', segCD);
    % D to A
    rt = (1-p_t).*ptD + p_t .* ptA;
    fld_t = subs(Av, [rx, ry, rz], [rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(3), p_t, 1)]; 
    integrand = dot(fld_t, dl);
    segDA = int(integrand, p_t, 0, 1);
    fprintf('DA segment: %.1f\n', segDA);
    fprintf('\t\t total flux = %.1f\n', segAB + segBC + segCD + segDA);
    
    %confirm by total differential
    % A to B rx[1,0] , ry = 1, rz = 0
    dl = [1, 0, 0];
    integrand = dot(Av, dl);
    integrand = subs(integrand, [ry, rz], [1, 0]);
    fluxAB = int(integrand, rx, 1, 0);
    fprintf('\nAB flux: %.1f\n', fluxAB);
    % B to C rx = 0 , ry = 1, rz[0, 2]
    dl = [0, 0, 1];
    integrand = dot(Av, dl);
    integrand = subs(integrand, [rx, ry], [0, 1]);
    fluxBC = int(integrand, rz, 1, 0);
    fprintf('BC flux: %.1f\n', fluxBC);
    % C to D rx = [0, 1] , ry = 1, rz = 2
    dl = [1, 0, 0];
    integrand = dot(Av, dl);
    integrand = subs(integrand, [ry, rz], [1, 2]);
    fluxCD = int(integrand, rx, 0, 1);
    fprintf('CD flux: %.1f\n', fluxCD);
    % D to A rx = 1 , ry = 1, rz[2, 0]
    dl = [0, 0, 1];
    integrand = dot(Av, dl);
    integrand = subs(integrand, [rx, ry], [1, 1]);
    fluxDA = int(integrand, rz, 2, 0);
    fprintf('DA flux: %.1f\n', fluxDA);
    fprintf('\t\ttotal flux = %.1f\n', fluxAB + fluxBC + fluxCD + fluxDA); 
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    % you could use the sum of surface integrals...but divergance theorem is much better
    Dv = [rz*rx^2, ry^3, ry*rz^2];  % rx[-1,1] , ry[0,4] , rz[1,3]
    divgD = ee.getDivgRec(Dv);
    verts = [ 1, 4, 1;
              1, 4, 3;
              1, 0, 1;
              1, 0, 3;
              -1, 0, 1;
              -1, 0, 3;
              -1, 4, 1;
              -1, 4, 3];
        %fun_graphVF(Dv, verts, 3);
        
    % the divergance can be used as the integrand... just take it on d{rx} d{ry} d{rz}
    flux = int(int(int(divgD, rx, -1, 1), ry, 0, 4), rz, 1, 3);
    fprintf('flux on this volume is: %.2f\n', flux);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    Av = [sr, -3, 5*sf];
    divgA = ee.getDivgSph(Av);
        %{
        ee.help_diff();
        display(divgA);           TRANSFORMS WORK
        test = [2, pi/4, pi/4];
        divgA = subs(divgA, [sr, st, sf], test);
        display(double(divgA));
        ee.feed(test, 'S');
        test = ee.pts(1,:);
        display(test);
        temp = ee.getU_SR();
        Avr = temp * transpose(Av);
        Avr = ee.simpleSRsub(Avr);
        Avr = simplify(transpose(Avr));	
        divg = ee.getDivgRec(Avr);
        divg = subs(divg, [rx, ry, rz], test);
        display(double(divg));
        %}

    % anyway, just take divergence sr[0,4] , st[0, pi/2] , sf[0, 3pi/2]
    % but keep in mind diff vol for spere: (sr)^2 sin(st) d{sr} d{st} d{sf}  SCALAR
    flux = int(int(int( (sr^2) * sin(st) * divgA, sr, 0, 4),...
            st, 0, sym(pi)/2),...
                sf, 0, 3*sym(pi)/2);
    fprintf('flux is: %s     , ...  %.3f\n', flux, double(flux));  % 484.584
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    intg1 = rx*ry; % over rx[0,1] , ry[0,1], rz[0,2]
    v1 = int(int(int(intg1, rx, 0, 1),...
            ry, 0, 1),...
                rz, 0, 2);
    fprintf('volume1 is: %.2f\n', v1);
    
    % for cyn volume, rememeber =  cr d{cr} d{cf} d{cz}  SCALAR
    intg2 = cr*cz;
    v2 = int(int(int(cr * intg2, cr, 1, 3),...
            cf, 0, sym(pi)),...
                cz, 0, 2);
    fprintf('volume2 is: %s   .... %.2f\n', v2, v2);
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    v1 = 6*rx*ry - 2*rx*rz + rz;
    v1grad = ee.getGradRec(v1);
    fprintf('grad of v1 : [ %s , %s, %s ]\n', v1grad);
    
    v2 = 10*cr*cos(cf) - cr*cz;
    v2grad = ee.getGradCyn(v2);
    fprintf('grad of v2 : [ %s , %s, %s ]\n', v2grad);
    
    v3 = (2*cos(sf))/sr;
    v3grad = ee.getGradSph(v3);
    fprintf('grad of v3 : [ %s , %s, %s ]\n', v3grad);
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    sclV = 10*rx*ry*rz - 2*rz*rx^2;
    ptV = [-1, 4, 3];
    gradV = ee.getGradRec(sclV);
    fprintf('the gradient of scalar V is [ %s , %s , %s ]\n', gradV);
    numV = subs(gradV, [rx, ry, rz], ptV);
    fprintf('value at point: [ %.2f , %.2f , %.2f ]\n', numV);
    
    sclU = 2*cr*sin(cf) + cr*cz;
    ptU = [2, sym(pi)/2, -1];
    gradU = ee.getGradCyn(sclU);
    fprintf('\nthe gradient of scalar U is [ %s , %s , %s ]\n', gradU);
    numU = subs(gradU, [cr, cf, cz], ptU);
    fprintf('value at point: [ %.2f , %.2f , %.2f ]\n', numU);
    
    sclW = (4/sr)*sin(st)*cos(sf);
    ptW = [ 1, sym(pi)/6, sym(pi)/2];
    gradW = ee.getGradSph(sclW);
    fprintf('\nthe gradient of scalar W is [ %s, , %s, %s ]\n', gradW);
    numW = subs(gradW, [sr, st, sf], ptW);
    fprintf('value at point: [ %s , %s, %s ]\n', numW);
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    rv = [rx, ry, rz]; % given a position vector of point, show nambla(r)^n = n(r^n-2) r
    syms n; % n is any integer
    
    rvMag = sqrt(sum( rv.^2));  % scalar (rx^2 + ry^2 + rz^2)^(1/2)
    rvMag_n = rvMag^n;%(rx^2 + ry^2 + rz^2)^(n/2) 
    grad_rvMag_n = ee.getGradRec(rvMag_n); % lhs
    
    rhs = n*(rvMag^(n-2)) .* rv;
    display(simplify(grad_rvMag_n - rhs)); % if their difference is 0, then they are the same  
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    sclT = rx^2 + ry^2 - rz;
    pt = [1, 1, 2];
    gradT = ee.getGradRec(sclT);
    path = subs(gradT, [rx, ry, rz], pt);
    fprintf(' take path of gradient... [ %.1f , %.1f , %.1f ]\n', path);
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    fam = rx-2*ry+rz; % family of planes, want unit normal 
    gradF = ee.getGradRec(fam);  % gradient is technically the normal vector  to any plane
    normU = gradF ./ norm(gradF);
    fprintf(' the unit normal vector to planes is   [ %.2f , %.2f , %.2f ]\n', normU);
    fprintf(" don't for get about\t\t       [ %.2f , %.2f , %.2f ]\n", -1 .* normU);
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    sclT = sr * sin(st) * cos(sf); % give sph scalar function
    ptT = [2, sym(pi)/30, sym(pi)/6];
    gradT = ee.getGradSph(sclT);
    val = subs(gradT, [sr, st, sf], ptT);
    fprintf('the gradient is [ %s , %s , %s ]\n', gradT);
    fprintf('value is [ %s , %s , %s ]\n', val);
    fprintf('.....  [ %.3f , %.3f , %.3f ]     mag and direction of max rate of change\n', val);
end


%------------------------------------------------------------------------------------------ #25
if sel == 25
    sclF = ry*rx^2 - 2*rx*ry^2 + rz^3; % given scalar function
    ptF = [2, 4, -3]; % given point
    dir = [1, 2, -1]; % given direction
    
    gradF = ee.getGradRec(sclF);              % find gradient
    nrm = subs(gradF, [rx, ry, rz], ptF);     % gradient makes normal vector
    dirU = dir ./ norm(dir);                  % unit vector out of direction
    dirDer = dot(nrm, dirU);  % directional derivative is projection of normal onto UNIT direction
    fprintf('the directional derivative at point is : %.2f\n', dirDer);
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    pln1 = rx + 2*ry + 2*rz - 5; % first plane
    grad1 = ee.getGradRec(pln1); % gradient of plane 1 provides normal vector to plane 1
    pln2 = rx + ry; % second plane
    grad2 = ee.getGradRec(pln2); % gradient of plane 2 provides normal vector to plane 2
    % now all you need is a cross product relationship to solve for angle
    quo = dot(grad1, grad2) / ( norm(grad1) * norm(grad2) );
    ang = acosd(quo);
    fprintf(' the angle between the 2 planes (by norma vector) : %.3f°  or %.3f rad\n',...
        ang, deg2rad(ang));
end


%------------------------------------------------------------------------------------------ #27
if sel == 27
    sclV = 4*rx*ry*exp(rz);
    ptV = [3, 1, -2];
    gradV = ee.getGradRec(sclV);
    gradVp = subs(gradV, [rx, ry, rz], ptV);
    fprintf('at point, max change on:  [ %.4f , %.4f , %.4f ]   mag = %.2f\n',...
        gradVp, norm(gradVp));
end
 


%------------------------------------------------------------------------------------------ #28
if sel == 28
    sclU = 3*rx*ry*rz;
    sclV = 5*ry*rx^2 + 2*ry*rz;
    lhs = ee.getGradRec( sclU * sclV);
    rhs = sclU .* ee.getGradRec(sclV) + sclV .* ee.getGradRec(sclU);
    display( simplify(lhs - rhs));  % match
end


%------------------------------------------------------------------------------------------ #29
if sel == 29
    
    Av = [ rx*ry, ry^2, -rx*rz];
    divgA = ee.getDivgRec(Av);
    fprintf('divg A = %s\n', divgA);
    
    Bv = [ cr*cz^2, cr*sin(cf)^2, 2*cr*cz*sin(cf)^2];
    divgB = ee.getDivgCyn(Bv);
    fprintf('\ndivg B = %s\n', divgB);
    
    Cv = [ sr, 0, sr*cos(st)^2];
    divgC = ee.getDivgSph(Cv);
    fprintf('\ndivg C = %s\n', divgC);
end



%------------------------------------------------------------------------------------------ #30
if sel == 30
    
    Av = [ ry*rx^2, rx, 2*ry*rz ];
    ptA = [-3, 4, 2];
    divgA = ee.getDivgRec(Av);
    valA = subs(divgA, [rx, ry, rz], ptA);
    fprintf('divg @ A =  %.3f\n', valA);
    
    Bv = [3*cr*sin(cf), -5*cz*cr^2, 8*cz*cos(cf)^2];
    ptB = [5, sym(pi)/6, 1];
    divgB = ee.getDivgCyn(Bv);
    valB = subs(divgB, [cr, cf, cz], ptB);
    fprintf('divg @ B =  %.3f\n', double(valB));
    
    Cv = [ cos(sf) * sr^2 , 0, 2*sr, ];
    ptC = [2, sym(pi)/3, sym(pi)/2];
    divgC = ee.getDivgSph(Cv);
    valC = subs(divgC, [sr, st, sf], ptC);
    fprintf('divg @ C =  %.3f\n', double(valC));
end


%------------------------------------------------------------------------------------------ #31
if sel == 31
    sclT = 50*sin( (sym(pi) * rx)/2 )  * cosh( sym(pi) * ry / 2);
    gradT = ee.getGradRec(sclT);
    fprintf('the gradient is  [ %s , %s , %s ] \n', gradT);
    divgT = ee.getDivgRec(gradT);
    fprintf('the divergence is  %s\n', divgT);
end


%------------------------------------------------------------------------------------------ #32
if sel == 32
    % trying to prove #330 vector to scalar divg relationship
    sclV = rx*ry*rz;
    Av = [ 2*rx, 3*ry, -4*rz];
    
    lhs = ee.getDivgRec( sclV .* Av );
    rhs = sclV * ee.getDivgRec(Av) + dot( Av, ee.getGradRec(sclV));
    fprintf('lhs: % s\n', simplify(lhs));
    fprintf('rhs: %s \n',simplify(rhs));
end
  

%------------------------------------------------------------------------------------------ #33
if sel == 33                            % be careful what side of nambla terms end up on
    rv = [rx, ry, rz];
    Tv = [2*rz*ry, rx*ry^2, ry*rz*rx^2];
    
    partA = ee.getDivgRec(rv) .* Tv;
    partA = simplify(partA);
    fprintf('A:  [ %s , %s, %s ]\n', partA); % [ 6*ry*rz , 3*rx*ry^2, 3*rx^2*ry*rz ]
    
    temp1 = rv(1) * diff(Tv, rx, 1);
    temp2 = rv(2) * diff(Tv, ry, 1);
    temp3 = rv(3) * diff(Tv, rz, 1);
    partB = temp1 + temp2 + temp3;
    fprintf('B:  [ %s , %s, %s ]\n', partB); % [ 4*ry*rz , 3*rx*ry^2, 4*rx^2*ry*rz ]
    
    temp1 = dot(rv, Tv);
    temp2 = ee.getDivgRec(rv);
    partC = temp1 * temp2;
    partC = simplify(partC);
    fprintf('C:  %s\n', partC);   % 3*rx*ry*(2*rz + rx*rz^2 + ry^2)
    
    temp = ee.getGradRec( norm(rv)^2 );
    partD = dot(temp, rv);
    fprintf('D:  %s    ...2r^2\n', simplify(partD));  %  2*rx^2 + 2*ry^2 + 2*rz^2    ...2r^2 
end
   

%------------------------------------------------------------------------------------------ #34
if sel == 34
    Av = [ 2*rx, -rz^2, 3*rx*ry];
    % want flux cr = 2, cf [0, pi/2] , cz[0,1]    implies convert to cyn
    
    temp = ee.getU_RC();
        %display(temp);
    temp = temp * transpose(Av);
        display(temp);
    temp = ee.simpleRCsub(temp);
    AvC = simplify(transpose(temp)); % the vector field is now in cylindrical
    fprintf('vector in cyn [ %s ,\n\t        %s , \n\t        %s ]\n', AvC);
    
    % use stokes...
    % A:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} > ...cr = 2
    
    %curlA = ee.getCurlCyn(AvC); no curl here... FLUX not CIRCULATION
    
    dS = [ cr, 0, 0];
    integrand = dot(AvC, dS);
    integrand = subs(integrand, cr, 2);
    integrand = simplify(integrand);
    fprintf('\nintegrating: % s    d{cf} d{cz}\n', integrand);
    flux = int(int(integrand, cf, 0, sym(pi)/2),...
            cz, 0, 1);
    fprintf('\nflux on surface:  %s    ...%.3f\n', flux, double(flux));
end


%------------------------------------------------------------------------------------------ #35
if sel == 35
    Dv = [ 2*cr*cz^2, 0, cr*cos(cf)^2];
    % cr [ 2, 5]  cf[0, 2pi], cz[-1,1]
    
    % top, bottom, inside, outside    sA: <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    %top, cr[2, 5], cf[0, 2pi], cz = 1
    dS = [ 0, 0, cr]; % cz is fixed
    intg = dot(Dv, dS);
    intg = subs(intg, cz, 1);
    topFlux = int(int(intg, cr, 2, 5),...
                    cf, 0, 2*sym(pi));
   fprintf('top flux is %s    ....%.3f\n', topFlux, topFlux);
   %bottom, cr[2, 5], cf[0, 2pi], cz = -1
   dS = [ 0, 0, -cr]; % cz is fixed, but orientate
   intg = dot(Dv, dS);
   intg = subs(intg, cz, -1);
   botFlux = int(int(intg, cr, 2, 5),...
                    cf, 0, 2*sym(pi));
   fprintf('bottom flux is %s    ....%.3f\n', botFlux, botFlux);
   %inside  cr = 2, cf[0, 2pi], cz[-1,1]
   dS = [-cr, 0, 0]; % cr fixed, but orientate
   intg = dot(Dv, dS);
   intg = subs(intg, cr, 2);
   isFlux = int(int(intg, cf, 0, 2*sym(pi)),...
                cz, -1, 1);
   fprintf('inside flux is %s    ....%.3f\n', isFlux, isFlux);
   %outside cr = 5, cf[0, 2pi], cz[-1, 1]
   dS = [cr, 0, 0]; % fixed on cr, orientation good
   intg = dot(Dv, dS);
   intg = subs(intg, cr, 5);
   osFlux = int(int(intg, cf, 0, 2*sym(pi)),...
                cz, -1, 1);
   fprintf('outside flux is %s    ....%.3f\n', osFlux, osFlux);
   fprintf('total flux: %s  ... %.3f\n',...
       topFlux + botFlux + isFlux + osFlux, double(topFlux + botFlux + isFlux + osFlux));
    
    % just use divg theorem  diff vol =  cr d{cr} d{cf} d{cz}  SCALAR
    divgD = ee.getDivgCyn(Dv);
    flux = int(int(int(cr * divgD, cr, 2, 5),...
                    cf, 0, 2*sym(pi)),...
                        cz, -1, 1);
    fprintf('\nflux is : %s    ... %.3f\n', flux, flux); % 56*pi  ... 175.929   much easier
end


%------------------------------------------------------------------------------------------ #36
if sel == 36
    Hv = [10*cos(st), 0, 0]; % flux on sr = 1, st[0, pi/2] , sf[0, 2pi]
    % sA:      <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    > 
    % vol:  (sr)^2 sin(st) d{sr} d{st} d{sf} 
    
    % top sr = 1, st[0, pi/2] , sf[0, 2pi]
    dS = [ sin(st) * sr^2, 0, 0]; % orientation is good
    intg = dot( Hv, dS);
    intg = subs(intg, sr, 1);
    topF = int(int(intg, st, 0, sym(pi)/2),...
                    sf, 0, 2*sym(pi));
    fprintf('top flux:  %s   ...  %.3f\n', topF, topF);
    % bot sr[0, 2], st = pi/2 , sf[0, 2pi]
    dS = [0, sr * sin(st), 0]; % orientation is good
    intg = dot(Hv, dS);
    intg = subs(intg, st, sym(pi)/2);
    botF = int(int(intg, sr, 0, 2),...
                    sf, 0, 2*sym(pi));
    fprintf('bot flux: %s    ... %.3f\n', botF, botF); % makes sense...cos(pi/2) = 0 on vector
    
    % divg way...
    divgH = ee.getDivgSph(Hv);
    %display(divgH);
    intg = sin(st) * sr^2 * divgH;
    %display(simplify(intg));
    flux = int(int(int(intg, sr, 0, 1),...
                st, 0, sym(pi)/2),...
                    sf, 0, 2*sym(pi));
    fprintf('\nflux is:  %s   ... %.3f\n', flux, flux);   % so much easier 10*pi   ... 31.416
end


%------------------------------------------------------------------------------------------ #37
if sel == 37
    Hv = [2*rx*ry, (rx^2) + (rz^2), 2*ry*rz]; % rx[0,1] , ry[1,2] , rz[-1,3]
    
    % easy way --> dig theorem
    divgH = ee.getDivgRec(Hv);
    flux = int(int(int(divgH, rx, 0, 1),...
                    ry, 1, 2),...
                        rz, -1, 3);
    fprintf('flux on box is: %.3f\n', flux);    % 24.000 it is
    
    % top, rz = 3
    dS = [0, 0, 1]; % orentation good
    intg = dot(Hv, dS);
    intg = subs(intg, rz, 3);
    f1 = int(int(intg, rx, 0, 1), ry, 1, 2);
    fprintf('\ntop flux: %.3f\n', f1);
    % bot, rz = -1
    dS = [0, 0, -1]; % orentation reversed
    intg = dot(Hv, dS);
    intg = subs(intg, rz, -1);
    f2 = int(int(intg, rx, 0, 1), ry, 1, 2);
    fprintf('bottom flux: %.3f\n', f2);
    % side1, x = 0
    dS = [-1, 0, 0]; % orentation reversed
    intg = dot(Hv, dS);
    intg = subs(intg, rx, 0);
    f3 = int(int(intg, ry, 1, 2), rz, -1, 3);
    fprintf('side1 flux: %.3f\n', f3);
    % side2, x = 1
    dS = [1, 0, 0]; % orentation good
    intg = dot(Hv, dS);
    intg = subs(intg, rx, 1);
    f4 = int(int(intg, ry, 1, 2), rz, -1, 3);
    fprintf('side2 flux: %.3f\n', f4);
    % side3, y = 1
    dS = [0, -1, 0]; % orentation reversed
    intg = dot(Hv, dS);
    intg = subs(intg, ry, 1);
    f5 = int(int(intg, rx, 0, 1), rz, -1, 3);
    fprintf('side4 flux: %.3f\n', f5);
    % side4, y = 2
    dS = [0, 1, 0]; % orentation good
    intg = dot(Hv, dS);
    intg = subs(intg, ry, 2);
    f6 = int(int(intg, rx, 0, 1), rz, -1, 3);
    fprintf('side4 flux: %.3f\n', f6);
    fprintf('total flux on volume: %.3f\n', f1+f2+f3+f4+f5+f6);   % lots of work
end


%------------------------------------------------------------------------------------------ #38
if sel == 38
    Hv = [ 4*cr^2, 0, -2*cz]; % cr = 10, cf[0,2pi] , cz[0,3]
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    % vol:  cr d{cr} d{cf} d{cz}
    % easy money way...
    divgH = ee.getDivgCyn(Hv);
    flux = int(int(int( cr * divgH, cr, 0, 10), cf, 0, 2*sym(pi)), cz, 0, 3);
    fprintf('flux on volume is:  %s   , ...  %.3f\n', flux, flux); % 23400*pi  ... 73513.268
    
    % top, cr[0,10], cf[0,2pi], cz = 3
    dS = [0, 0, cr]; % orientation good
    intg = dot(Hv, dS);
    intg = subs(intg, cz, 3);
    topF = int(int(intg, cr, 0, 10), cf, 0, 2*sym(pi));
    fprintf('\ntop flux:   %s   ... %.3f\n', topF, topF);
    % bot, cr[0,10], cf[0,2pi], cz = 0
    dS = [0, 0, -cr]; % orientation reversed
    intg = dot(Hv, dS);
    intg = subs(intg, cz, 0);
    botF = int(int(intg, cr, 0, 10), cf, 0, 2*sym(pi));
    fprintf('bot flux:   %s   ... %.3f\n', botF, botF);
    % side flux, cr = 10 , cf[0, 2pi], cz[0,3]
    dS = [cr, 0, 0]; % oreientation good
    intg = dot(Hv, dS);
    intg = subs(intg, cr, 10);
    sidF = int(int(intg, cf, 0, 2*sym(pi)), cz, 0, 3);
    fprintf('side flux:  %s  ... %.3f\n', sidF, sidF);
    fprintf('adding surfaces, total flux is:  %s  ... %.3f\n',...
        topF + botF + sidF, double(topF + botF + sidF));  
end


%------------------------------------------------------------------------------------------ #39
if sel == 39
    Av = [rx^2, ry^2, rz^2]; % cr = 1 , cz[2, 4] ...transform implied
    % you'll need to take this rec Av and make it cyn, find divg
    % then take it cr[0,1] , cf[0,2pi], cz[2,4]
    AvC = ee.transVecRC(Av);
    fprintf('transformed vector: [ %s , %s , %s ]\n', AvC);
    divgA = ee.getDivgCyn(AvC);
    divgA = simplify(divgA);
    fprintf('...and integrate this divergence: %s\n', divgA);
    flux = int(int(int(cr*divgA, cr, 0, 1), cf, 0, 2*sym(pi)), cz, 2, 4);
    fprintf('flux on volume is:   %s    ....  %.3f\n', flux, flux);
    
    %the savage way...
    dvg = ee.getDivgRec(Av);
    flx = int(int(int(dvg, rz, 2, 4),...
                ry, -sqrt(1-rx^2), sqrt(1-rx^2)),...
                    rx, -1, 1);
    fprintf(' \n with no param change flux:    %s    ....  %.3f\n', flx, flx);
end


%------------------------------------------------------------------------------------------ #40
if sel == 40
    Av = [ sr^2, sr * sin(st) * cos(sf), 0]; % sr[0,3] , st[0, pi/2] , sf[0, pi/2]
    dvg = ee.getDivgSph(Av);
    % vol:  (sr)^2 sin(st) d{sr} d{st} d{sf}
    intg = dvg * sin(st) * sr^2;
    temp1 = int(intg, sr, 0, 3);
    temp2 = int(temp1, st, 0, sym(pi)/2);
    flux = int(temp2, sf, 0, sym(pi)/2);
    fprintf('flux is:  %s   ...  %.3f\n', flux, flux);
    
    % the verificaiton...
    % sA:      <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    > 
    % top, sr = 3 , st[0, pi/2] , sf[0, pi/2]
    dS = [sin(st) * sr^2, 0, 0]; % orientation good
    intg = dot(Av, dS);
    intg = subs(intg, sr, 3);
    topF = int(int(intg, st, 0, sym(pi)/2), sf, 0, sym(pi)/2);
    fprintf('\ntop flux:   %s   ...  %.3f\n', topF, topF);
    % bottom flux, sr[0,3] , st = pi/2 , sf[0, pi/2
    dS = [0, sr * sin(st), 0]; % orientation is good
    intg = dot(Av, dS);
    intg = subs(intg, st, sym(pi)/2);
    botF = int(int(intg, sr, 0, 3), sf, 0, sym(pi)/2);
    fprintf('bot flux:   %s   ...  %.3f\n', botF, botF);
    % face 1, sr[0, 3] , st[0,pi/2] , sf = 0
    dS = [0, 0, -sr]; % reversed
    intg = dot(Av, dS);
    intg = subs(intg, sf, 0);
    f1f = int(int(intg, sr, 0, 3), st, 0, sym(pi)/2);
    fprintf('face1 flux:   %s   ...  %.3f\n', f1f, f1f);
    
    % face 2, sr[0, 3] , st[0,pi/2] , sf = pi/2
    dS = [0, 0, sr]; % good
    intg = dot(Av, dS);
    intg = subs(intg, sf, sym(pi)/2);
    f2f = int(int(intg, sr, 0, 3), st, 0, sym(pi)/2);
    fprintf('face2 flux:   %s   ...  %.3f\n', f2f, f2f);
    fprintf('add surface fluxes  =  %s   ...  %.3f\n',...
        topF + botF + f1f + f2f, double(topF + botF + f1f + f2f));  
    % match (81*pi)/2 + 9   ...  136.235
end


%------------------------------------------------------------------------------------------ #41
if sel == 41
    Fv = [sin(cf)*cr^2, cz*cos(cf), cr*cz]; % cr[2,3] , cf[0, 2pi] , cz[0,5]
    FvR = ee.transVecCR(Fv);
        %fun_graphVF(FvR, 0, 2);
    % find flux of big obj and subtrac little one?
    %vol:  cr d{cr} d{cf} d{cz}
    divgF = ee.getDivgCyn(Fv);
    smallF = int(int(int( cr*divgF, cr, 0, 2), cf, 0, 2*sym(pi)), cz, 0, 5);
    fprintf('interior flux: %s   ...  %.3f\n', smallF, smallF);
    bigF = int(int(int( cr*divgF, cr, 0, 3), cf, 0, 2*sym(pi)), cz, 0, 5);
    fprintf('exterior flux: %s   ...  %.3f\n', bigF, bigF);
    netF = bigF - smallF;
    fprintf('net flux: %s   ...  %.3f\n', netF, netF);
   
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} > 
    % top cr[2, 3] , cf[0, 2pi] , cz = 5
    dS = [0, 0, cr]; % orientation good
    intg = dot(Fv, dS);
    intg = subs(intg, cz, 5);
    topF = int(int(intg, cr, 2, 3), cf, 0, 2*sym(pi));
    fprintf('\ntop flux:  %s   ... %.3f \n', topF, topF);
    % bottom cr[2, 3] , cf[0, 2pi] , cz = 0
    dS = [0, 0, -cr]; % orientation reversed
    intg = dot(Fv, dS);
    intg = subs(intg, cz, 0);
    botF = int(int(intg, cr, 2, 3), cf, 0, 2*sym(pi));
    fprintf('bot flux:  %s   ... %.3f \n', botF, botF);
    % inside cr = 2 , cf[0, 2pi] , cz [0, 5]
    dS = [-cr, 0, 0]; % reversed
    intg = dot(Fv, dS);
    intg = subs(intg, cr, 2);
    inF = int(int(intg, cf, 0, 2*sym(pi)), cz, 0, 5);
    fprintf('inside flux =  %s  .. %.3f\n', inF, inF);
    % outside cr = 3 , cf[0, 2pi] , cz [0, 5]
    dS = [cr, 0, 0]; % good orientation
    intg = dot(Fv, dS);
    intg = subs(intg, cr, 3);
    outF = int(int(intg, cf, 0, 2*sym(pi)), cz, 0, 5);
    fprintf('outside flux =  %s  .. %.3f\n', outF, outF);
    fprintf('total: %s   ...  %.3f\n',...
        topF + botF + inF + outF, double(topF + botF + inF + outF));
    % (190*pi)/3   ...  198.968    subtract big flux from little one and call it good
end


%------------------------------------------------------------------------------------------ #42
if sel == 42
    Av = [rx*ry, ry^2, -rx*rz];
    curlA = ee.getCurlRec(Av);
    fprintf('curl of A: [ %s , %s , %s ]\n', curlA);
    
    Bv = [cr*cz^2, cr*sin(cf)^2, 2*cr*cz*sin(cf)^2];
    curlB = simplify(ee.getCurlCyn(Bv));
    fprintf('curl of B: [ %s , %s , %s ]\n', curlB);
    
    Cv = [sr, 0, sr*cos(st)^2];
    curlC = simplify(ee.getCurlSph(Cv));
    fprintf('curl of B: [ %s , %s , %s ]\n', curlC);
end


%------------------------------------------------------------------------------------------ #43
if sel == 43
    % looking at curl and divg(curl)
    Av = [ ry*rx^2, rz*ry^2, -2*rx*rz];
    curlA = ee.getCurlRec(Av);
    dcurlA = ee.getDivgRec(curlA);
    fprintf('curl is [ %s, %s, %s ] ...  divg(curl) = %s\n', curlA, dcurlA);
    
    Av = [ cz*cr^2, cr^3, 3*cr*cz^2];
    curlA = ee.getCurlCyn(Av);
    dcurlA = ee.getDivgCyn(curlA);
    fprintf('\ncurl is [ %s, %s, %s ] ...  divg(curl) = %s\n', curlA, dcurlA);
    
    Av = [sin(sf)/sr^2, -cos(sf)/sr^2,0];
    curlA = simplify(ee.getCurlSph(Av));
    dcurlA = ee.getDivgSph(curlA);
    fprintf('\ncurl is [ %s, %s, %s ] ...  divg(curl) = %s\n', curlA, dcurlA);
end


%------------------------------------------------------------------------------------------ #44
if sel == 44
    Hv = [ cr*sin(cf), cr*cos(cf), -cr]; 
    curlH = simplify(ee.getCurlCyn(Hv));
    fprintf('the curl is [ %s , %s, %s ]\n', curlH);
    curlCurlH = simplify(ee.getCurlCyn(curlH));
    fprintf('the curl of the curl is [ %s , %s, %s ]\n', curlCurlH);
end


%------------------------------------------------------------------------------------------ #45
if sel == 45
    denom = ( rx^2 + ry^2 + rz^2)^(3/2);
    Gv = [ rx/denom, ry/denom, rz/denom]; % the gravity vector field
        %fun_graphVF(Gv, 0, 3);
    curlG = ee.getCurlRec(Gv);
    fprintf('curl of gravity vector is [ %s, %s, %s ]\n', curlG);
end


%------------------------------------------------------------------------------------------ #46
if sel == 46
    Fv = [ry*rx^2, -ry, 0];
    ptA = [0, 0, 0];
    ptB = [1, 1, 0];
    ptC = [2, 0, 0];
        %fun_graphVF(Fv, [ptA;ptB;ptC], 0);
    % by stokes
    curlF = ee.getCurlRec(Fv);
    dS = [ 0, 0, -1]; % they are trying to be tricking, boundary oriented downward
    intg = dot(curlF, dS);
    intg = subs(intg, rz, 0);
    circ1 = int(int(intg, ry, 0, rx), rx, 0, 1); % going to split it up into 2 regions
    circ2 = int(int(intg, ry, 0, 2-rx), rx, 1, 2);
    circ = circ1 + circ2;
    fprintf('stokes says:  %s    ...  %.3f\n', circ, circ);   % 7/6  ...  1.167
    
    % line integral way   segment 1
    rt = (1-p_t) .* ptA + p_t * ptB;
    Fv_t = subs(Fv, [rx, ry, rz], [ rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(2), p_t, 1) ];
    intg = dot(Fv_t, dl);
    seg1 = int(intg, p_t, 0, 1);
    fprintf('\nseg1:   %s   ...   %.3f\n', seg1, seg1);
    % line integral way  segment 2
    rt = (1-p_t) .* ptB + p_t * ptC;
    Fv_t = subs(Fv, [rx, ry, rz], [ rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(2), p_t, 1) ];
    intg = dot(Fv_t, dl);
    seg2 = int(intg, p_t, 0, 1);
    fprintf('seg2:   %s   ...   %.3f\n', seg2, seg2);
    % line integral way  segment 3
    rt = (1-p_t) .* ptC + p_t * ptA;
    Fv_t = subs(Fv, [rx, ry, rz], [ rt(1), rt(2), rt(3)]);
    dl = [ diff(rt(1), p_t, 1), diff(rt(2), p_t, 1), diff(rt(2), p_t, 1) ];
    intg = dot(Fv_t, dl);
    seg3 = int(intg, p_t, 0, 1);
    fprintf('seg3:   %s   ...   %.3f\n', seg3, seg3);
    seg = seg1 + seg2 + seg3;
    fprintf(' total:  %s  ...  %.3f\n   stokes confirmed\n', seg, seg);
end


%------------------------------------------------------------------------------------------ #47
if sel == 47
    % there is a right way and a wrong way to do this. cyn is the right way ...rec is savage
    Fv = [ry*rx^2, -ry, 0]; % going to put in cyn and go cr[1,2] , cf[0,pi/2], cz=0
    Fv_cyn = ee.transVecRC(Fv);
    fprintf('vf trans to cyn: [ %s , %s , %s ]\n', Fv_cyn);  % could even do sph if you want
    curlF = ee.getCurlCyn(Fv_cyn);
    curlF = simplify(curlF);
    fprintf('integrate this curl [stokes] :  [ %s , %s, %s ]\n', curlF);
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    dS = [ 0, 0, -cr]; % need orientation adjusted
    intg = dot(curlF, dS);
    intg = subs(intg, cz, 0); % just a formality
    circ = int(int(intg, cr, 1, 2), cf, 0, sym(pi)/2);
    fprintf('stokes says: %s   ...  %.3f\n', circ, circ);
    
    % length:  <  d{cr}           ,  cr d{cf}     ,  d{cz}          >
    %seg1 cr[1,2] , cf = pi/2 , cz = 0
    dl = [1, 0, 0];
    intg = dot(Fv_cyn, dl);
    intg = subs(intg, [cf, cz], [sym(pi)/2, 0]);
    seg1 = int(intg, cr, 1, 2);
    fprintf('\nseg1 = %s ...  %.3f\n', seg1, seg1);
    
    %seg2 cr=2 , cf[pi/2, 0] , cz = 0
    dl = [0, cr, 0];  % opposes, but don't change...limits of integration
    intg = dot(Fv_cyn, dl);
    intg = subs(intg, [cr, cz], [2, 0]);
    seg2 = int(intg, cf, sym(pi)/2, 0);
    fprintf('seg2 = %s ...  %.3f\n', seg2, seg2);
    
    %seg3 cr[2,1] , cf = 0 , cz = 0
    dl = [1, 0, 0];
    intg = dot(Fv_cyn, dl);
    intg = subs(intg, [cf, cz], [0, 0]);
    seg3 = int(intg, cr, 2, 1);
    fprintf('seg3 = %s ...  %.3f\n', seg3, seg3);
    
    %seg4 cr=1 , cf[0, pi/2] , cz = 0
    dl = [0, cr, 0];
    intg = dot(Fv_cyn, dl);
    intg = subs(intg, [cr, cz], [1, 0]);
    seg4 = int(intg, cf, 0, sym(pi)/2);
    fprintf('seg4 = %s ...  %.3f\n', seg4, seg4);
    circ = seg1+seg2+seg3+seg4;
    fprintf('circ is:  %s   ... %.3f\n', circ, circ);  % match  (15*pi)/16   ... 2.945
end


%------------------------------------------------------------------------------------------ #48
if sel == 48
    Av = [cr * sin(cf), cr^2, 0]; % using same figure, but no transform needed now
    % take it cr[1,2] , cf[0, pi/2], cz=0
    curlA = ee.getCurlCyn(Av);
    % sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    % by stokes...
    dS = [0, 0, -cr]; % need to orientate
    intg = dot(curlA, dS);
    intg = subs(intg, cz, 0); % formality
    circ = int(int(intg, cr, 1, 2), cf, 0, sym(pi)/2);
    fprintf('stokes says: %s      ...  %.3f\n', circ, circ);   %  3/2 - (7*pi)/2  ...  -9.496
    
    % length:  <  d{cr}           ,  cr d{cf}     ,  d{cz}          >
    %seg1 cr[1,2] , cf = pi/2 , cz = 0
    dl = [1, 0, 0];
    intg = dot(Av, dl);
    intg = subs(intg, [cf, cz], [sym(pi)/2, 0]);
    seg1 = int(intg, cr, 1, 2);
    fprintf('\nseg1:   %s   ...  %.3f\n', seg1, seg1);
    %seg2 cr = 2 , cf = [ pi/2 , 0 ] , cz = 0
    dl = [0, cr, 0];
    intg = dot(Av, dl);
    intg = subs(intg, [cr, cz], [2, 0]);
    seg2 = int(intg, cf, sym(pi)/2, 0);
    fprintf('seg2:   %s   ...  %.3f\n', seg2, seg2)
    %seg3 cr[2,1] , cf = 0 , cz = 0
    dl = [1, 0, 0];
    intg = dot(Av, dl);
    intg = subs(intg, [cf, cz], [0, 0]);
    seg3 = int(intg, cr, 2, 1);
    fprintf('seg3:   %s   ...  %.3f\n', seg3, seg3);
    %seg4 cr = 1 , cf = [ 0, pi/2 ] , cz = 0
    dl = [0, cr, 0];
    intg = dot(Av, dl);
    intg = subs(intg, [cr, cz], [1, 0]);
    seg4 = int(intg, cf, 0, sym(pi)/2);
    fprintf('seg4:   %s   ...  %.3f\n', seg4, seg4);
    circ = seg1+seg2+seg3+seg4;
    fprintf('total is    %s  ...  %.3f\n', circ, circ);
end


%------------------------------------------------------------------------------------------ #49
if sel == 49
    Fv = [ 2*cr*cz, 3*cz*sin(cf), -4*cr*cos(cf) ]; % verify on cr[0,2] , cf[0,pi/4] , cz = 1
    %  sA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    curlF = ee.getCurlCyn(Fv);
    dS = [0, 0, cr];
    intg = dot(curlF, dS);
    intg = subs(intg, cz, 1);
    circ = int(int(intg, cr, 0, 2), cf, 0, sym(pi)/4);
    fprintf('stokes says circ is:  %s  ... %.3f\n', circ, circ);%6 - 3*2^(1/2)  ... 1.757
    
    % length:  <  d{cr}           ,  cr d{cf}     ,  d{cz}          >   to verify
    %seg1 cr[0,2] , cf = 0 , cz = 1
    dl = [1, 0, 0];
    intg = dot(Fv, dl);
    intg = subs(intg, [cf, cz], [0, 1]);
    seg1 = int(intg, cr, 0, 2);
    fprintf('\nseg1: %s ... %.3f\n', seg1, seg1);
    
    %seg2 cr=2, cf[0, pi/4] , cz = 1
    dl = [0, cr, 0];
    intg = dot(Fv, dl);
    intg = subs(intg, [cr, cz], [2, 1]);
    seg2 = int(intg, cf, 0, sym(pi)/4);
    fprintf('seg2: %s ... %.3f\n', seg2, seg2);
    
    %seg3 cr[2,0] , cf = pi/4 , cz = 1
    dl = [1, 0, 0];
    intg = dot(Fv, dl);
    intg = subs(intg, [cf, cz], [sym(pi)/4, 1]);
    seg3 = int(intg, cr, 2, 0);
    fprintf('seg3: %s ... %.3f\n', seg3, seg3);
    circ = seg1+seg2+seg3;
    fprintf('\ttotal:  %s  ... %.3f\n', circ, circ); % have to keep same direction as surface
end


%------------------------------------------------------------------------------------------ #50
if sel == 50
    Av = [4*exp(-ry)*rx^2, -8*rx*exp(-ry), 0];
    temp1 = ee.getDivgRec(Av); % the divergence of Av
    temp2 = ee.getGradRec(temp1); % gradient of divg
    sol = ee.getCurlRec(temp2); % and take the curl
    fprintf(' sol:  [ %s , %s , %s ]\n', simplify(sol) );
end


%------------------------------------------------------------------------------------------ #51
if sel == 51
    sclV = sin(st)*cos(sf)/sr;
    gradV = simplify( ee.getGradSph(sclV) );
    fprintf('grad V :   [ %s , %s, %s ]\n', gradV);
    curlGradV = simplify( ee.getCurlSph(gradV) );
    fprintf(' curl(gradV) :  [ %s, %s, %s ]\n', curlGradV);
    divgGradV = simplify( ee.getDivgSph(gradV) );
    fprintf('divg(gradV) :  %s\n', divgGradV);
    check = ee.getLaplacianSphS(sclV);
    fprintf('               %s\n', simplify(check)); % should be same as scalar Laplacian
end
    

%------------------------------------------------------------------------------------------ #52
if sel == 52
    mult = sqrt( rx^2 + ry^2 + rz^2 ) / sqrt( rx^2 + ry^2);
    Qv = [ (rx - ry ) * mult, ( rx + ry ) * mult, 0 ]; % going to want a sph conversion
        %fun_graphVF(Qv, [0, 0, sqrt(3)], sqrt(3));
    QvS = ee.transVecRS(Qv);
    fprintf('transformed vector: [ %s ,\n\t\t      %s,\n\t\t      %s ]\n', simplify(QvS)); % looks ok
    
    %{
        length:  <  d{sr}                        ,   sr d{st}              , sr sin(st) d{sf}  >  
        sA:      <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    >  
        vol:  (sr)^2 sin(st) d{sr} d{st} d{sf}  
    %}
    
    % the line integral sr = 2, st = pi/6, sf[0, 2pi] by regular line...
    dl = [ 0, 0, sr * sin(st) ];
    intg = dot(QvS, dl);
    intg = subs(intg, [sr, st], [2, sym(pi)/6]);
    circ = int(intg, sf, 0, 2*sym(pi));
    fprintf('\nA by line: %s .... %.3f\n', circ, circ);  % 4*pi .... 12.566  &&&
    cynQ = ee.transVecRC(Qv); % cyn isn't bad... cr[0,1], cf[0,2pi], cz = sqrt(3)
    crrlQ = simplify(ee.getCurlCyn(cynQ)); %  ... rec integral crashes
    dS = [ 0, 0, cr];
    intg = dot(crrlQ, dS);
    intg = subs(intg, cz, sqrt(3));
    circ = int(int(intg, cr, 0, 1), cf, 0, 2*sym(pi));
    fprintf(' stokes on cyn  %s ... %.3f\n', circ, circ);
    
    % s2 is -4pi because it is reversed...still on same line integral
    
    % a normal sphere would be rx^2 + ry^2 + rz^2 = sr^2  ....this sphere is:
    %  rx^2 + ry^2 + (rz - sqrt(3))^2 = 1
    % implies regular rx, ry, rz transform, but not sr = radius equation
    sphEqnRec = rx^2 + ry^2 + (rz - sqrt(3))^2;
    sphEqnSph = simplify(ee.simpleRSsub(sphEqnRec));
    sphEqnSph = subs(sphEqnSph, sr, 1);  %  cos(st) = sqrt(3) / 2 
    
    % s1 surface flux  sr = 1 , st = [0, pi/6] , sf[0,2pi]
    QvC = ee.transVecRC(Qv);
    dS = [(sr)^2 * sin(st), 0, 0];
    intg = dot(QvS, dS)
    intg = subs(intg, sr, 2);
    s1f = int(int(intg, st, 0, sym(pi)/2), sf, 0, 2*sym(pi));
    fprintf('flux on s1:   %s   ... %.3f \n', s1f);
    
    
    
    % s2 surface flux  cr[0, 2] , cf[0, 2pi], cz[0, sqrt(3) ] 
    
    
    % volume int --> divg is easy  sr[2, 3] , st[0, pi/6] , sf[0,2pi]
    divgQs = simplify(ee.getDivgSph(QvS));
    intg = sin(st) * (cr^2) * divgQs;
    flux_sph = int(int(int(sr, sqrt(3), 2+sqrt(2)), st, 0, sym(pi)/6), sf, 0, 2*sym(pi));
    flux_cone = int(int(int(sr, 0, sqrt(3)), st, 0, sym(pi)/6), sf, 0, 2*sym(pi));
    fprintf(' %s   ... %.3f\n', simplify(flux_sph + flux_cone) , flux_sph + flux_cone);
end

%------------------------------------------------------------------------------------------ #53
if sel == 53
    Hv = [ 2*rx*rz, 5*rx*ry*rz, 8*(ry + rz)];
    divgH = ee.getDivgRec(Hv);
    curlH = ee.getCurlRec(Hv);
    fprintf('divg = %s\n', simplify(divgH));
    fprintf('curl = [ %s , %s , %s ]\n', simplify(curlH));
end


%------------------------------------------------------------------------------------------ #54
if sel == 54
    Bv = [sr^2, 4*sr*cos(2*st), 0];
    divgB = ee.getDivgSph(Bv);
    curlB = ee.getCurlSph(Bv);
    fprintf('divg is: %s\n', simplify(divgB));
    fprintf('\ncurl is: [ %s , %s , %s ]\n', simplify(curlB));
end


%------------------------------------------------------------------------------------------ #55
if sel == 55
    recV = rx*ry*rz;
    gradV = ee.getGradRec(recV);
    temp = recV .* gradV;
    lhs = ee.getDivgRec(temp);
    rhs = recV * ee.getLaplacianRecS(recV) + ( norm(gradV) )^2;
    fprintf(' rhs:  %s\n', simplify(rhs));
    fprintf(' lhs:  %s\n', simplify(lhs)); % it worked...they are equal
    
    Av = [ rx, ry, rz];
    lhs = ee.getCurlRec( recV .* Av);
    rhs = recV .* ee.getCurlRec(Av) + cross(gradV, Av);
    fprintf('\n rhs:  [ %s , %s , %s ]\n', simplify(rhs));
    fprintf(' lhs:  [ %s , %s , %s ]\n', simplify(lhs)); % it worked...they are equal
end


%------------------------------------------------------------------------------------------ #56
if sel == 56
    Bv = [ry*rx^2, ry+2*rx^2, rz-ry]; % given
    dvg = ee.getDivgRec(Bv);
    fprintf('divg: %s\n', dvg);
    crl = ee.getCurlRec(Bv);
    fprintf('curl: [ %s , %s , %s ]\n', crl);
    gradDivg = ee.getGradRec(dvg);
    fprintf('grad of divg: [ %s, %s, %s ]\n', gradDivg);
    crlCrl = ee.getCurlRec(crl);
    fprintf('curl of curl: [ %s , %s, %s ]\n', crlCrl);
    temp = curl(Bv)
end


%------------------------------------------------------------------------------------------ #57
if sel == 57
    v1 = rx^3 + ry^3 +rz^3;
    v1L = ee.getLaplacianRecS(v1);
    fprintf(' %s \n', simplify(v1L));
    v2 = sin(2*cf)*cr*cz^2;
    v2L = ee.getLaplacianCynS(v2);
    fprintf(' %s \n', simplify(v2L));
    v3 = (sr^2) * ( 1 + cos(st)*sin(sf) );
    v3L = ee.getLaplacianSphS(v3);
    %fprintf(' %s \n', simplify(v3L));
    pretty(simplify(v3L));
end


%------------------------------------------------------------------------------------------ #58
if sel == 58
    
    U = (rx^3)*(ry^2)*exp(rx*rz);
    ptU = [1, -1, 1];
    lapU = ee.getLaplacianRecS(U);
    valU = subs(lapU, [rx, ry, rz], ptU);
    fprintf('laplacian:  %s\n', simplify(lapU));
    fprintf('value:  %s  ___  %.3f\n', valU, valU);
    
    V = (cr^2)*cz*(cos(cf) + sin(cf));
    ptV = [5, sym(pi)/6, -2];
    lapV = ee.getLaplacianCynS(V);
    valV = subs(lapV, [cr, cf, cz], ptV);
    fprintf('\nlaplacian:  %s\n', lapV);
    fprintf('value:  %s  ___  %.3f\n', valV, double(valV));
    
    W = exp(-sr)*sin(st)*cos(sf);
    ptW = [1, sym(pi)/3, sym(pi)/6];
    lapW = ee.getLaplacianSphS(W);
    valW = subs(lapW, [sr, st, sf], ptW);
    %fprintf('\nlaplacian:  %s\n', simplify(lapW));
    pretty(simplify(lapW));
    fprintf('value:  %s  ___  %.3f\n', valW, valW);
end
    

%------------------------------------------------------------------------------------------ #59
if sel == 59
    rv = [rx, ry, rz];
    r = norm(rv);
    lhs = ee.getGradRec( log(r) );
    rhs = rv ./ r^2;
    fprintf('lhs1:  %s \n', simplify(lhs(1)));
    fprintf('rhs1:  %s \n', simplify(rhs(1)));
    fprintf('\n\tlhs2:  %s \n', simplify(lhs(2)));
    fprintf('\trhs2:  %s \n', simplify(rhs(2)));
    fprintf('\n\t\tlhs3:  %s \n', simplify(lhs(3)));
    fprintf('\t\trhs3:  %s \n', simplify(rhs(3)));
    
    dump = input('any key:', 's');
    clc;
    
    lhs = ee.getLaplacianRecS( log(r) );
    rhs = 1 / (r^2);
    lhsV = subs(lhs, [rx, ry, rz], [ 1, 2, 3]);
    rhsV = subs(rhs, [rx, ry, rz], [ 1, 2, 3]);
    fprintf('lhs1:  %s \n', simplify(lhs(1)));
    fprintf('rhs1:  %s \n', simplify(rhs(1)));
    fprintf(' lhs: %s  ,  rhs: %s\n', lhsV, rhsV);
end


%------------------------------------------------------------------------------------------ #60
if sel == 60
    U = rx*(ry^2)*(rz^3);
    gradU = ee.getGradRec(U);
    testU = ee.getDivgRec(gradU); % divergence of gradient is same as laplacian
    lapU = ee.getLaplacianRecS(U);
    fprintf('grad:  [ %s , %s , %s ]\n', gradU);
    fprintf('lap: %s\n', lapU);
    
    scl = sin(cf)/cr;
    grad = ee.getGradCyn(scl);
    test = ee.getDivgCyn(grad);
    lap = ee.getLaplacianCynS(scl);
    fprintf('\ngradient: [ %s , %s , %s ] \n', grad);
    fprintf('laplacian: %s\n', lap);
    
    scl = (sr^2)*sin(st)*cos(sf);
    grad = ee.getGradSph(scl);
    test = ee.getDivgSph(grad);
    lap = ee.getLaplacianSphS(scl);
    temp = simplify(lap - test);
    fprintf('\ngradient: [ %s , %s , %s ] \n', grad);
    fprintf('laplacian: %s\n', simplify(lap));
end


%------------------------------------------------------------------------------------------ #61
if sel == 61
    sclr = (cr^2)*cz*cos(cf);
    grad = ee.getGradCyn(sclr);
    lap = ee.getLaplacianCynS(sclr);
    fprintf('grad: [ %s , %s , %s ]\n', grad);
    fprintf('lap: %s \n', lap);
end


%------------------------------------------------------------------------------------------ #62
if sel == 62
    sclr = 5*cos(sf)/(sr^2);
    grad = ee.getGradSph(sclr)
    lap = ee.getLaplacianSphS(sclr)
    curlGrad = ee.getCurlSph(grad)
end

%------------------------------------------------------------------------------------------ #63
if sel == 63
    sclr = 4*rx*ry*rz^2 + 10*ry*rz;
    lap = ee.getLaplacianRecS(sclr)
    comp = ee.getDivgRec( ee.getGradRec(sclr))
end

%------------------------------------------------------------------------------------------ #64
if sel == 64
    Av = [ 2*cr*sin(cf), 4*cr*cos(cf), cr + cr*cz^2];
    vecLap = ee.getLaplacianCynV(Av);
    pretty(simplify(vecLap));
end

%------------------------------------------------------------------------------------------ #65
if sel == 65
    Av = [rx*rz, rz^2, ry*rz];
    temp0 = ee.getDivgRec(Av);
    lhs = simplify( ee.getGradRec( temp0));
    pretty(lhs);
    temp1 = ee.getDivgRec(Av);
    temp2 = ee.getGradRec(temp1);
    temp3 = ee.getLaplacianRecV(Av);
    rhs = simplify( temp2 - temp3);
    pretty(rhs);
end


%------------------------------------------------------------------------------------------ #66
if sel == 66
    Av = [rx, ry, rz];
        %fun_graphVF(Av, [0,0,2], 2);
    divg = ee.getDivgRec(Av);
    crl = ee.getCurlRec(Av);
    fprintf('\ndivg = %s   ,    divg not 0 --> NOT solenoidal\n', simplify(divg));
    fprintf('curl = [ %s , %s , %s ]  ,   curl = [ 0, 0, 0 ] --> IS irrational/potential/conservative\n', simplify(crl));
    
    Bv = [2*cr*cos(cf), -4*cr*sin(cf), 3];
        temp = ee.transVecCR(Bv);
        %fun_graphVF(temp, [0,0,2], 2);
    divg = ee.getDivgCyn(Bv);
    crl = ee.getCurlCyn(Bv);
    fprintf('\ndivg = %s   ,    divg = 0 --> IS solenoidal\n', simplify(divg));
    fprintf('curl = [ %s , %s , %s ]  ,   curl not [ 0, 0, 0 ] --> NOT irrational/potential/conservative\n', simplify(crl));
    
    Cv = [ sin(st), 0, sr*sin(st)];
        temp = ee.transVecSR(Cv);
        %fun_graphVF(temp, [0,0,2], 2);
    divg = ee.getDivgSph(Cv);
    crl = ee.getCurlSph(Cv);
    fprintf('\ndivg = %s   ,    divg not 0 --> NOT solenoidal\n', simplify(divg));
    fprintf('curl = [ %s , %s , %s ]  ,   curl not [ 0, 0, 0 ] --> NOT irrational/potential/conservative\n', simplify(crl));
end


%------------------------------------------------------------------------------------------ #67
if sel == 67
    Gv = [ 16*rx*ry-rz, 8*rx^2, -rx];
        %fun_graphVF(Gv, [0,0,2],2);
    dvg = ee.getDivgRec(Gv);
    crl = ee.getCurlRec(Gv);
    fprintf('divg: %s       , not 0 , NOT solenoidal\n', dvg);
    fprintf('curl:   [ %s , %s , %s ] , if curl = [ 0, 0, 0 ] --> IS  irrational/potential/conservative\n', crl);
    
    % use divg theorem for net flux on cube  rx[0, 1], ry[0, 1] , rz[0, 1]
    intg = ee.getDivgRec(Gv);
    tempZ = int(intg, rz, 0, 1);
    tempY = int(tempZ, ry, 0, 1);
    flux = int(tempY, rx, 0,1);
    fprintf('\nnext flux unit cube:  %s  ___ %.3f\n', flux, flux );
    
    % flux on square rx[0,1] , ry[0,1], z = 0 
    % stokes or add lines  keep anti-clock...normal orientation
    ptA = [0,0,0];
    ptB = [1,0,0];
    ptC = [1,1,0];
    ptD = [0,1,0];
    temp = ee.getCurlRec(Gv);
    dS = [0, 0, 1];
    intg = dot( temp, dS);
    tempX = int(intg, rx, 0, 1);
    circ = int(tempX, ry, 0, 1);
    fprintf('circulation on unit square:  %s    _____  %.3f\n', circ, circ); 
end


%------------------------------------------------------------------------------------------ #68
if sel == 68
    E = [ 1/(2*sym(pi)*cr) , 0 , 0];
    Er = ee.transVecCR(E);
        %fun_graphVF(Er, [0,0,2],2);
    crl = ee.getCurlCyn(E);
    dvg = ee.getDivgCyn(E);
    fprintf('divg = %s      solenoidal\n', dvg);
    fprintf('curl = [ %s , %s , %s ]   conservative\n', crl);
    f = potential([10/rx^2,0,0],[rx,ry,rz])
end
    






%{
    
    Q = ee.transVecRS(Qv);
    divgQ = ee.getDivgSph(Q);
    intg = simplify(divgQ * sin(st) * sr^2);
    flux = int(int(int(intg, sr, sqrt(3), 2), st, 0, sym(pi)/6), sf, 0, 2*sym(pi));
    fprintf('\nvol flux:  %.3f\n', flux);
    
    sphV = [sr * (cos(sf) - sin(sf)), sr * (cos(sf) + sin(sf)), 0];
    temp = ee.getU_RS();
    sphV = simplify( temp * transpose(sphV));
    sphV = transpose(sphV);
    divV = (sr^2) * sin(st) * ee.getDivgSph(sphV);
    volF = int(int(int(divV, sr, 2, sqrt(3)), st, 0, sym(pi)/6), sf, 0, 2*sym(pi));
    fprintf('vol flux :   %.3f\n',  volF);
    display(double(volF+flux));
%}
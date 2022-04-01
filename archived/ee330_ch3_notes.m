%{
    chapter 3 is the last of the preliminary work ... should be good if you can handle it
    
    #1   ex3.1
    #2   pp3.1
    #3   ex3.2
    #4   pp3.2
    #5   ex3.3     and good derivation of namblas
    #6   pp3.3     gradient from scalar function
    #7   ex3.4     directional derivative
    #8   pp3.4     directional derivative
    #9   ex3.5     a line and an ellipsoid...holy shit
    #10  pp3.5     two functions, two normals, find angle
    #11  ex3.6     divergence
    #12  pp3.6     divg and eval @ pt
    #13  ex3.7     divg to find out-flux
    #14  pp3.7     divg to find out-flux, still on a cylinder
    #15  ex3.8     curl
    #16  pp3.8     curl at point
    #17  ex3.9     stokes on surface
    #18  pp3.9     stokes
    #19  ex3.11    laplacians on scalars
    #20  pp3.11    laplacians on scalars
    #21  vector fields demonstartion     4 main types  defined by divergence and curl
    #22  pp3.12     prove conservation   in two ways  ...curl(A) = [0,0,0]


%}
clc;
close all;
clearvars;


                       sel = 22;  % CHANGE CHANGE CHANGE


ee = cls_EE330_helper();
global rx; global ry; global rz; % rectangular params  rx, ry, rz        
global cr; global cf; global cz; % cylindrical params  cr, cf, cz
global sr; global st; global sf; % spherical params    sr, st, sf
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); syms cf; assume(cf,'real'); syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); syms st; assume(st,'real'); syms sf; assume(sf, 'real');
global par_t; syms par_t; assume(par_t, 'real'); % for the space curve

%------------------------------------------------------------------------------------------ #1
if sel == 1
    % given an obvious cylinder
    ptAr = [5, 0, 0]; % given point in rec
    ee.feed(ptAr, 'R');
    ptAc = ee.pts(2,:);
    ptBr = [0, 5, 0];
    ee.feed(ptBr, 'R');
    ptBc = ee.pts(2,:);
    ptCr = [0, 5, 10];
    ee.feed(ptCr, 'R');
    ptCc = ee.pts(2,:);
    ptDr = [5, 0, 10];
    ee.feed(ptDr, 'R');
    ptDc = ee.pts(2,:);
    
    myRad = 5; myCent = [0, 0]; myHt = [0, 10]; myPts = [ptAr; ptBr; ptCr; ptDr; 0, 0, 0];
    fun_graph_cyn('x',myRad, myCent, myHt, myPts)
    %(axis,rad, cent, hgt, pts)
    
    % length BC can be done in rec
    lenBC = sqrt(sum( (ptBr - ptCr).^2 ) );
    fprintf('from B to C is: %.3f\n', lenBC);
    % this is equiv to integrating [0, 10] w/r cz  or rz
    
    %[func] = @(x) x=1;
    %lenBC = integral(func, 0, 10);
    %display(lenBC);
    
    lenBC = int(1, cz, 0, 10); % litterarlly just integrate IAW displacement vector
    fprintf('from B to C is: %.3f\n', lenBC); 
    
    %lengthCD   notice dl = cr d(cf)
    lenCD = sqrt(sum( (ptCr - ptDr).^2 ));  
    fprintf('\nfrom C to D is: %.3f\n', lenCD);
    lenCD = int(5, cf, 0, pi/2);  % cr is fixed at 5,   w/r cf  [0, pi/2]
    fprintf('from C to D is: %.3f\n', lenCD);  % notice difference...it is not a straigt line
    
    % surface area ABCD has a normal of acr    diff is cr d(cf) d(cz)   
    % you can see cr will stay fixed at 5, cf [0, pi/2] , cz [0, 10]
    intCF = int(5, cf, 0, pi/2);
    intCZ = int(intCF, cz, 0, 10);
    saABCD = intCZ;
    fprintf('\nsa ABCD, normal to acr: %.3f\n', saABCD);
    
    % surface area ABO has a normal -acz   diff is  cr d(cr) d(cf)
    % no need to consider -cz for volume, just need surface area
    % cr [0, 5] and cf [0, pi/2] on this surface
    intCR = int(cr, cr, 0, 5);
    saAOB = int(intCR, cf, 0, pi/2);
    fprintf('\nsa AOB, anti-normal to acz: %.3f\n', saAOB);
    
    % surface AOFD is anti acf, but that is ok here diff = d(cr) d(cz)
    % with cr [0, 5] and cz[0, 10]
    intCR = int(1, cr, 0, 5);
    saAODF = int(intCR, cz, 0, 10);
    fprintf('\nsa AODF, anti-normal to acf: %.3f\n', saAODF); % makes sense... 5 * 10 = 50
    
    % and the volume.... diff =  cr d(cr) d(cf) d(cz)   ranges [0, 5], [0, pi/2], [0, 10]
    intCR = int(cr, cr, 0, 5);
    intCF = int(intCR, cf, 0, pi/2);
    vol = int(intCF, cz, 0, 10);
    fprintf('\nvolume: %.3f\n', vol);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    % sr [ 3, 5]    st [pi/3, pi/2]   sf [pi/4, pi/3]
    ee.diff_ping();
    
    % arc length DH    going D [3, pi/2, pi/4]  to H [3, pi/2, pi/3]  
    lenDH = int(3 * sin(pi/2), sf, pi/4, pi/3);  % differential on asf
    fprintf('lenDH: %.4f\n', lenDH);
    
    % arc length FG   going G [ 5, pi/2, pi/3 ]   to F [ 5, pi/3, pi/3 ]
    % only a differential length change on  ast ...[pi/3 , pi/2]
    % diff = sr  d(st)
    lenFG = int(5, st, pi/3, pi/2);
    fprintf('lenFG: %.4f\n',lenFG);
    
    % surf AEHD is in the asr direction  (anti) , but just want surface area
    % diff = (sr)^2  sin(st)  d(st) d(sf)
    intST = int(3^2 * sin(st), st, pi/3, pi/2);
    saAEHD = int(intST, sf, pi/4, pi/3);
    fprintf('AEHD area : %.4f\n', saAEHD);
    
    % surf ABCD is taking the asf direction 
    % diff =  sr d{sr} d{st}      sr[3, 5] , st[pi/3 , pi/2]
    intSR = int(sr, sr, 3, 5);
    saABCD = int(intSR, st, pi/3, pi/2);
    fprintf('ABCD area: %.4f\n', saABCD);
    
    % volume is easy... (sr)^2 sin(st) d{sr} d{st} d{sf}  on given interval ranges
    intSR = int(sr^2, sr, 3, 5);
    intST = int(sin(st), st, pi/3, pi/2); % using that seperate and multiply technique
    intSF = int(1, sf, pi/4, pi/3);
    vol = intSR * intST * intSF;
    fprintf('vol: %.4f\n', vol); 
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % calculating circulation /  line integral along a closed path
    ee.diff_ping();
    Fv = [rx^2, -rx*rz, -ry^2]; % given
    
    ptList = [1, 0, 0;
              0, 0, 0;
              0, 1, 0;
              1, 1, 1];
    fun_graph_path3(ptList, Fv); % (pts, Vf)
    %fun_graphVF(Fv,0);
    %fun_graph_sph(3, [0,0,0],9);  %(rad, cent, pts)
    %fun_graphVF([rz*cos(ry),rz^2,4*rx], [0,0,0]); % (field, pts)
    %fun_graph_norm(ry^2 - rx^2, [0,1,0]); % (funIn, pts)
    
    % part 1 is a segment from [1, 0, 0] to [0, 0, 0]   here, y=z=0
    % always take dl in the +a direction
    dl1 = [1, 0, 0];
    integrand1 = dot(Fv, dl1); 
    seg1 = int(integrand1, rx, 1, 0); % the path will produce a negative anwer
    fprintf('first integral = %.3f\n', seg1);
    
    % on segment 2, you can see dl = [0, 1, 0] dy   but Fv = [0,0,-ry^2]
    % dp = 0, so integral is 0
    seg2 = 0;
    fprintf('second integral = %.3f\n', seg2);
    
    % on segment 3, dl = [1, 0, 0] dx  + [0, 0, 1]dz   ,  Fv = [ rx^2, -rx*rz , -1 ]
    % ry does not change...consistent with definition
    integrand3x = rx^2; %  dot ( [ rx^2, -rx*rz , -1 ] ,  [1, 0, 0] ) dx
    integrand3z = -1;  % dot (  [ rx^2, -rx*rz , -1 ] ,  [0, 0, 1] ) dz
    seg3x = int(integrand3x, rx, 0, 1);
    seg3z = int(integrand3z, rz, 0, 1);
    seg3 = seg3x + seg3z;
    fprintf('third integral = %.3f\n', seg3);
    
    % segment 4 does not have a change in x ... but it is backwards opposed to ary and arz
    % if x is fixed at 1, FV = [ 1, -rz, -ry^2 ]
    % dl: [0, -1, 0] dy  + [0, 0, -1]dz  =>  [ 0, dy, dz ]
    % integrate  0 dx , rz dy , ry^2 dz
    % you could have also kept standard ary arz unit vectors and seen 1 -> 0 (not 0 ->1)
    % also, ry = rx on this path, so you can combine
    % int [1, 0] -rz - ry^2  = [0, 1] rz + ry^2  and sub ry = rz
    integrand4y = -ry;
    integrand4z = -ry^2;
    seg4 = int(integrand4y+integrand4z, ry, 1, 0);
    fprintf('fourth integral = %.3f\n', seg4);
    
    fprintf('final answer: %.3f\n', seg1+seg2+seg3+seg4);
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    ee.diff_ping();
    % a wedge in cyn   cr[0, 2] , cf[0,pi/3] , cz = 0
    Av = [cr * cos(cf), 0, cz*sin(cf)]; % given      want circulation
    
    %seg1, d{cz} is out of the picture, so is d{cf} it is all d{cp}
    % if cz = 0 and cf = 0, Av = [ cr, 0, 0 ]
    % dl = [ d{cr}, 0, 0 ] for this path and it goes [0, 2]
    % ig1 = dot( Av, dl) = cr d{cr}
    seg1 = int( cr, cr, 0, 2);
    fprintf('first segment = %.3f\n', seg1);
    
    %seg2, we know cz = 0 the entire time, and the radius does not increase cr = 2
    % dl = [0, cr d{cf}, 0]  where cf[0, pi/3]
    % Av = [ 2*cos(cf) , 0, 0]
    % ig2 = dot(Av, dl) = 0
    seg2 = 0;
    fprintf('second segment = %.3f\n', seg2);
    
    %seg3 , cz still = 0, cr[0,2], and cf does not change, = pi/3
    % dl = [ d{cr}, 0, 0]
    % Av = [cr (1/2) , 0, 0]
    %ig3 = (1/2) cr d{cr} but take it [2,0] to keep signs right
    seg3 = int( (1/2)*cr, cr, 2, 0);
    fprintf('third segment = %.3f\n', seg3);
    
    fprintf('the total circulation is %.3f\n', seg1+seg2+seg3);
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    %ee.help_namb(); % still not seeing the nambla -> cyn, nambla -> sph
    
    scalV = exp(-rz)*sin(2*rx)*cosh(ry);
    gradV = gradient(scalV, [rx, ry, rz]);
    gradV = ee.getNambRec() .* transpose(gradV);
    fprintf('\n');
    fprintf('grad V = %s\n', gradV);
    dumpV = ee.getGradRec(scalV);
    
    scalU = cr^2 * cz * cos(2*cf);
    gradU = gradient(scalU, [cr, cf, cz]);
    gradU = ee.getNambRec() .* transpose(gradU);
    fprintf('\n');
    fprintf('grad U = %s\n', gradU);
    dumpU = ee.getGradCyn(scalU);
    
    scalW = 10*sr*cos(sf)*(sin(st))^2;
    gradW = gradient(scalW, [sr, st, sf]);
    gradW = ee.getNambRec() .* transpose(gradW);
    fprintf('\n');
    fprintf('grad W = %s\n', gradW);
    dumpW = ee.getGradSph(scalW);
end


%------------------------------------------------------------------------------------------ #6
if sel == 6                             % do these on paper...simplification not always gauruntee
    scalU = rx^2 * ry + rx*ry*rz;
    gradU = ee.getGradRec(scalU);
    fprintf('gradient: %s\n', gradU);
    
    scalV = cr*cz*sin(cf) + cz^2 * (cos(cf))^2 + cr^2;
    gradV = ee.getGradCyn(scalV);
    fprintf('gradient: %s\n', simplify(gradV));
    
    scalF = cos(st) * sin(sf) * log(sr) + sf*sr^2;
    gradF = ee.getGradSph(scalF);
    fprintf('gradient: %s\n', gradF);
    fprintf('\ngradient: %s\n', simplify(gradF));
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    scalW = rx^2 * ry^2 + rx*ry*rz;   % given 3var function
    pt = [2, -1, 0];  % given point
    vec = [3, 4, 12]; % given vector...find directional derivative
    
    % find gradW
    gradW = ee.getGradRec(scalW);
    % multiply this by the differential length [1, 1, 1] ...no change
    % now obtain the numeric gradient at the point
    gradWnum = subs(gradW, [rx, ry, rz], pt);
    display(gradWnum);
    
    % the directional derivative w/r to vec is simply the dot product of 
    %      numeric gradient 'gradWnum' and the unit vector of 'vec'
    vecU = vec ./ norm(vec);
    dWdl = dot(gradWnum, vecU);
    display(dWdl);
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    scalF = rx*ry + ry*rz + rx*rz;  % given 3 var scalar function
    tgt_pt = [1, 2, 3]; % given point to evaluate at
    d_pt = [3, 4, 4]; % given point to take direction to
    vec = d_pt - tgt_pt;  % implies this vector between points
    
    gradF = ee.getGradRec(scalF);
    gradFnum = subs(gradF, [rx, ry, rz], tgt_pt);  % evaluate gradient at target point given
    display(gradFnum);
    vecU = vec ./ norm(vec); % make a unit vector in desired direction
    dFdl = dot(gradFnum, vecU);
    display(dFdl); % boom   
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % an ellipsoid rx^2 + ry^2 + 2*rz^2 = 10  ...what angle does x=y=2z intersect?
    % or   (1/2) x = (1/2) y = z       
    spaceC = [2*param_t, 2*param_t, param_t]; % saying x = (1/2) t = y,  z = t
    ellipsd = (rx^2)/10 + (ry^2)/10 + (rz^2/5);%   x^2 / 10 + y^2/10 + z^2/5 = 1
    scalE = (rx^2) + (ry^2) + 2*(rz^2) - 10; % taking ellipsoid as scalar = 0  f(x,y,z) = 0
    tRange = [-5, 5];
    
    centr = [0,0,0];
    sAx = [sqrt(10), sqrt(10), sqrt(5)]; % implied axes of ellipsoid
    %fun_graph_spaceC(spaceC, tRange);  %(params,range)
    %fun_graphEllipsoid(centr,sAx); %(cent,semAx)
    
    % by the equation, if line is [2t, 2t, t]
    %     x = y = 2t  ,  z = t         so they meet at   (2t)^2 + (2t)^2 + 2*t^2 = 10
    % that means t^2 = 1,  so t = +/- 1
    pt1 = subs(spaceC, param_t, -1);
    display(pt1);
    pt2 = subs(spaceC, param_t, 1);  % using this one
    display(pt2);
    vecU = pt2 ./ norm(pt2);
    
    
    temp = ee.getGradRec(scalE); % find gradient of ellipsoid ... must used implied scalar
    gradE = subs(temp, [rx, ry, rz], pt2);  % and find it at t = -1
    display(gradE);
    % the angle is just a dp calculation of these 2 constant vectors
    normU = gradE ./ norm(gradE); % also the anti normal
    quo = dot(gradE, pt2) / abs(dot(gradE, pt2));
    
    display(double(quo));
    
    % make gradient a unit normal vector
    %gradEu = gradE ./ norm(gradE);
    %display(double(gradEu));
    % solve 
    %quo = dot(gradE, vecU) / norm(gradE);
    %ang = acos(quo);
    %display(double(rad2deg(quo)));
    %quo = dot(gradE, vec) / norm(gradE);  % no need for abs(vecU)....it is 1
    %ang = acosd(quo);
    %display(double(quo));
    % the line is point t = -1 to t = 1 through ellipsoid
    %vec = pt2 - pt1;  % so make a vector from t = -1 to t = 1
    %display(vec);
    %vecU = vec ./ norm(vec);  % and make the unit vector
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    fun1 = rx^2 * ry + rz - 3;
    fun2 = rx * log(rz) - ry^2 + 4;
    pt = [-1, 2, 1]; % it checks out as an intersection point
    temp = ee.getGradRec(fun1);
    grad1 = subs(temp, [rx, ry, rz], pt);
    temp = ee.getGradRec(fun2);
    grad2 = subs(temp, [rx, ry, rz], pt);
    fprintf(' at point, grad1 = [ %.1f, %.1f, %.1f ]\n', grad1);
    fprintf(' at point, grad2 = [ %.1f, %.1f, %.1f ]\n', grad2);
    
    % just a dp 
    quo = dot(grad1, grad2) / ( norm(grad1) * norm(grad2) );
    ang = acosd(quo);
    fprintf('the angle between two normals is %.3f°\n', 180-ang); % smalles angle
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    %ee.help_divg();
    Pv = [rx^2*ry*rz, 0, rx*rz];   %  rec VF
    divgP = ee.getDivgRec(Pv);
    fprintf('divg of Pv: %s\n', divgP);
    
    Qv = [cr*sin(cf), cr^2 * cz, cz * cos(cf)];  % cyn VF
    divgQ = ee.getDivgCyn(Qv);
    fprintf('\ndivg of Qv: %s\n', divgQ);
    
    Tv = [ (1/sr^2)*cos(st), sr*sin(st)*cos(sf), cos(st)]; % sph VF
    divgT = ee.getDivgSph(Tv);
    fprintf('\ndivg of Tv: %s\n', divgT);
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    Av = [ry*rz, 4*rx*ry, ry]; % rec VF
    ptA = [1, -2, 3];
    divgA = ee.getDivgRec(Av);
    fprintf('divg A: %s\n', divgA);
    divgApt = subs(divgA, [rx, ry, rz], ptA);
    fprintf('divg A at point: %.3f\n', divgApt);
    
    Bv = [cr*cz*sin(cf), 3*cr*(cz^2)*cos(cf), 0]; % cyn VF
    ptB = [5, pi/2, 1];
    divgB = ee.getDivgCyn(Bv);
    fprintf('\ndivg B: %s\n', divgB);
    divgBpt = subs(divgB, [cr, cf, cz], ptB);
    fprintf('divg B at point: %.3f\n', divgBpt);
    
    Cv = [2*sr*cos(st)*cos(sf),0, sqrt(sr)]; % sph VF
    ptC = [1, pi/6, pi/3];
    divgC = ee.getDivgSph(Cv);
    fprintf('\ndivg C: %s\n', divgC);
    divgCpt = subs(divgC, [sr, st, sf], ptC);
    fprintf('divg C at point: %.3f\n', divgCpt);
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    Gv = [cr, 0, 1] .* (10*exp(-2*cz));   % given VF on cylinder [0, 1], [0, 2pi], [0, 1]
    % there is no cf component
    % good thing it is a closed surface... just get divg and integrate with proper dv
    divgG = ee.getDivgCyn(Gv);
    fprintf('divg = %.2f      ...not going to be any flux\n', divgG);
          % no outward flux .... they cancel
    
    % if you had to do this without divergence theorem....top, bottom, and sides needed
    
    % on top, cr[0, 1] , cf[0, 2pi] , z = 1
    % surface area is in cz direction  (only place outward flux can occur
    % and diff is  <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    dS_top = [0, 0, cr];
    topIntg = dot(Gv, dS_top); 
    topIntg = subs(topIntg, cz, 1); % don't forget cz = 1
    fprintf('\ntop integrand cz = 1, int(cr, 0, 1) --> int(cf, 0, 2pi):  %s\n', topIntg);
    temp = int(topIntg, cr, 0, 1);
    topInt = int(temp, cf, 0, 2*pi);
    fprintf('the top surface flux is: %.3f\n', topInt);
    
    % on bottom, z = 0 , but you have to use -acz   ...flux opposes direction of acz
    dS_bot = [0, 0, -cr];
    botIntg = dot(Gv, dS_bot);
    botIntg = subs(botIntg, cz, 0);
    fprintf('\nbottom integrand cz = 0, int(cr, 0, 1) --> int(cf, 0, 2pi):  %s\n', botIntg);
    temp = int(botIntg, cr, 0, 1);
    botInt = int(temp, cf, 0, 2*pi);
    fprintf('the bottom surface flux is: %.3f\n', botInt);
    
    % on sides, cr = 1 , cf[0, 2pi] , cz[0, 1]  
    % ... only concerned about [1, 0, 0] direction   , what is emminating from sides in cr direction
    dS_sides = [ cr, 0, 0 ]; % with differential
    sidesIntg = dot(Gv, dS_sides);
    sidesIntg = subs(sidesIntg, cr, 1);
    fprintf('\nsides integrand cr = 0, int(cz, 0, 1) --> int(cf, 0, 2pi):  %s\n', botIntg);
    temp = int(sidesIntg, cz, 0, 1);
    sidesInt = int(temp, cf, 0, 2*pi);
    fprintf('the sides surface flux is: %.3f\n', sidesInt);
    
    % now just add them all up to confirm
    fprintf('\nthe net surface flux is: %.3f\n', topInt+botInt+sidesInt);
    
    temp = ee.getU_CR();
    Gv_rec = temp * transpose(Gv);
    Gv_rec = transpose(Gv_rec);
    Gv_rec = ee.simpleCRsub(Gv_rec);
    testv = [rz*rx^2, ry*ry*rz, rx+2];
    %fun_graphVF(testv, 0);
    display(Gv_rec);
    fun_graph_cyn('z',1, [0,0], [0,1],0);
end

%------------------------------------------------------------------------------------------ #14
if sel == 14
    % another cylinder cr[0, 4] , cf[0, 2pi] , cz[0,1]
    %fun_graph_cyn('z',4, [0,0], [0,1],0);
    Dv = [(cr^2)*(cos(cf))^2, cz*sin(cf), 0]; % cyn vecf
    divgD = ee.getDivgCyn(Dv);
    fprintf('divg = %s  ...it is not 0, there will be some flux\n', divgD);
    % do a quick check....diff volume is cr d{cr} d{cf} d{cz}
    integrand = divgD * cr;
    temp1 = int(integrand, cr, 0, 4);
    temp2 = int(temp1, cf, 0, 2*sym(pi));  % very slick
    fluxD = int(temp2, cz, 0, 1);
    fprintf('out-flux should be: %s\n', fluxD);
    
    % to do the top, cr[0, 4] , cf[0, 2pi] , z = 1
    dS = [0, 0, cr]; % outward flux can only occur in this direction
    integrand = dot(Dv, dS); % for surface area and flux across it
    integrand = subs(integrand, cz, 1); % fix z = 1
    integrand = int(integrand, cr, 0, 4);
    integral_top = int(integrand, cf, 0, 2*sym(pi));
    fprintf('\ntop outward-flux: %s\n', integral_top);
    
    % to do the bottom, cr[0, 4] , cf[0, 2pi] , z = 0
    dS = [0, 0, -cr]; % outward flux can only occur in this direction...anti
    integrand = dot(Dv, dS); % for surface area and flux across it
    integrand = subs(integrand, cz, 0); % fix z = 0
    integrand = int(integrand, cr, 0, 4);
    integral_bottom = int(integrand, cf, 0, 2*sym(pi));
    fprintf('bottom outward-flux: %s\n', integral_bottom);
    
    % to do the sides, cr = 4, cf[0, 2pi], cz[0, 1]
    dS = [cr, 0, 0]; % outward flux must emminate from direction of radius for side
    integrand = dot(Dv, dS);
    integrand = subs(integrand, cr, 4);
    integrand = int(integrand, cf, 0, 2*sym(pi));
    integral_sides = int(integrand, cz, 0, 1);
    fprintf('sides outward-flux: %s\n', integral_sides);
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    %ee.help_curl();
    Pv = [rx^2*ry*rz, 0, rx*rz];   %  rec VF
    curlP = ee.getCurlRec(Pv);
    fprintf('\ncurl of P[1] = %s\n', curlP(1));
    fprintf('curl of P[2] = %s\n', curlP(2));
    fprintf('curl of P[3] = %s\n', curlP(3));
    
    Qv = [cr*sin(cf), cr^2 * cz, cz * cos(cf)];  % cyn VF
    curlQ = ee.getCurlCyn(Qv);
    fprintf('\ncurl of Q[1] = %s\n', curlQ(1));
    fprintf('curl of Q[2] = %s\n', curlQ(2));
    fprintf('curl of Q[3] = %s\n', simplify(curlQ(3))); % be careful with shit that isn't simplified
    
    Tv = [ (1/sr^2)*cos(st), sr*sin(st)*cos(sf), cos(st)]; % sph VF
    curlT = ee.getCurlSph(Tv);
    fprintf('\ncurl of T[1] = %s\n', (curlT(1)));
    fprintf('curl of T[2] = %s\n', simplify(curlT(2)));
    fprintf('curl of T[3] = %s\n', simplify(curlT(3)));  % it works, but you have to simplify further
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    Av = [ry*rz, 4*rx*ry, ry]; % rec VF
    ptA = [1, -2, 3];
    curlA = ee.getCurlRec(Av);
    fprintf('\ncurl of A[1] = %s\n', simplify(curlA(1)));
    fprintf('curl of A[2] = %s\n', simplify(curlA(2)));
    fprintf('curl of A[3] = %s\n', simplify(curlA(3)));
    curlA_pt = subs(curlA, [rx, ry, rz], ptA);
    fprintf('curl at point: [ %.2f, %.2f, %.2f ]\n', curlA_pt);
    
    Bv = [cr*cz*sin(cf), 3*cr*(cz^2)*cos(cf), 0]; % cyn VF
    ptB = [5, pi/2, 1];
    curlB = ee.getCurlCyn(Bv);
    fprintf('\ncurl of B[1] = %s\n', simplify(curlB(1)));
    fprintf('curl of B[2] = %s\n', simplify(curlB(2)));
    fprintf('curl of B[3] = %s\n', simplify(curlB(3)));
    curlB_pt = subs(curlB, [cr, cf, cz], ptB);
    fprintf('curl at point: [ %.2f, %.2f, %.2f ]\n', curlB_pt);
    
    Cv = [2*sr*cos(st)*cos(sf),0, sqrt(sr)]; % sph VF
    ptC = [1, pi/6, pi/3];
    curlC = ee.getCurlSph(Cv);
    fprintf('\ncurl of C[1] = %s\n', simplify(curlC(1)));
    fprintf('curl of C[2] = %s\n', simplify(curlC(2)));
    fprintf('curl of C[3] = %s\n', simplify(curlC(3)));
    curlC_pt = subs(curlC, [sr, st, sf], ptC);
    fprintf('curl at point: [ %.2f, %.2f, %.2f ]\n', curlC_pt);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    V = [cr*cos(cf), sin(cf), 0]; % given cyn vec
    % want circulation on path    cr[2, 5]    cf[pi/6, pi/3]   cz[0,0]
    % lets check with stokes before getting into the line integrals
    curlV = ee.getCurlCyn(V);
    fprintf('\ncurl of V[1] = %s\n', simplify(curlV(1)));
    fprintf('curl of V[2] = %s\n', simplify(curlV(2)));
    fprintf('curl of V[3] = %s\n', simplify(curlV(3)));
    % this will also be the integrand...
        %keeping in mind area diff cyn = <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
            % curl is in acz direction   ...the normal to the surface
    dS = [0, 0, cr];    % and cr[2, 5]  cf[pi/6, pi/3]
    temp = dot(curlV, dS);
    integrand = subs(temp, cz, 0); % make sure and cz is set to 0
    temp = int(integrand, cr, 2, 5);
    flux = int(temp, cf, pi/6, pi/3);
    fprintf('the flux should be: %.3f\n', flux);
    
    %doing all the silly line integrals
    % keep in mind diff length is <  d{cr}           ,  cr d{cf}     ,  d{cz}          >
    % AB is cr = 2, cf[pi/3, pi/6] , cz = 0
    dl = [0, cr, 0]; % cf is the only direction that is changing
    integrand = dot(V, dl);
    integrand = subs(integrand, [cr, cz], [2, 0]);
    ABflux = int(integrand, cf, pi/3, pi/6);    % your limits of integration handle direction change
    fprintf('\nflux on side AB : %.3f\n', ABflux); %   everything is good as long as left hand on surf
    
    % BC is cr = [2, 5], cf = pi/6 , cz = 0
    dl = [1, 0, 0]; % cr is the only direction that is changing
    integrand = dot(V, dl);
    integrand = subs(integrand, [cf, cz], [pi/6, 0]);
    BCflux = int(integrand, cr, 2, 5);
    fprintf('flux on side BC : %.3f\n', BCflux);
    
    % CD is cr = 5 , cf = [ pi/6, pi/3 ] , cz = 0
    dl = [0, cr, 0]; % cf is the only direction that is changing
    integrand = dot(V, dl);
    integrand = subs(integrand, [cr, cz], [5, 0]);
    CDflux = int(integrand, cf, pi/6, pi/3);
    fprintf('flux on side CD : %.3f\n', CDflux);
    
    % DA is cr = [5, 2] , cf =  pi/3  , cz = 0
    dl = [1, 0, 0]; % cr is the only direction that is changing
    integrand = dot(V, dl);
    integrand = subs(integrand, [cf, cz], [pi/3, 0]);
    DAflux = int(integrand, cr, 5, 2);
    fprintf('flux on side DA : %.3f\n', DAflux);
    
    fprintf('the line integrals sum to: %.3f\n', ABflux + BCflux + CDflux + DAflux);
    fprintf('...perfect match, but stokes is better\n');
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    Av = [cr * cos(cf), 0, cz*sin(cf)]; % given
    % ... cr[0, 2] , cf[0, pi/3] , cz = 0
    curlA = ee.getCurlCyn(Av);
    fprintf('\ncurl of A[1] = %s\n', simplify(curlA(1)));
    fprintf('curl of A[2] = %s\n', simplify(curlA(2)));
    fprintf('curl of A[3] = %s\n', simplify(curlA(3)));
    % not that the gurl is ready, notice curl is only in z direction
    % diff cyn = <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >
    dS = [0, 0, cr]; % or [ cr, 1, cr]...won't matter
    integrand = dot(curlA, dS);
    integrand = subs(integrand, cz, 0); % just in case
    temp = int(integrand, cr, 0, 2);
    flux = int(temp, cf, 0, pi/3);
    fprintf('flux accross surface is : %.2f\n', flux);
end
  

%------------------------------------------------------------------------------------------ #19
if sel == 19
    %ee.help_laplac()
    V = exp(-rz)*sin(2*rx)*cosh(ry); % scalar field
    laplacV = ee.getLaplacianRecS(V);
    fprintf('namb2 of V : %s\n', simplify(laplacV));
    
    U = (cr^2) * cz * cos(2*cf);
    laplacU = ee.getLaplacianCynS(U);
    fprintf('\nnamb2 of U : %s\n', simplify(laplacU));
    
    W = 10 * sr * cos(sf) * (sin(st))^2;
    laplacW = ee.getLaplacianSphS(W);
    fprintf('\nnamb2 of U : %s\n', simplify(laplacW));
    pretty(simplify(laplacW));
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    U = (rx^2)*ry + rx*ry*rz; % rec scalar field
    laplacU = ee.getLaplacianRecS(U);
    fprintf('namb2 of V : %s\n', simplify(laplacU));
    
    V = cr * cz * sin(cf) + (cz^2)*(cos(cf))^2 + cr^2;  % cyn scalar field
    laplacV = ee.getLaplacianCynS(V);
    fprintf('\nnamb2 of U : %s\n', simplify(laplacV));
    
    F = cos(st)*sin(sf)*log(sr) + sf*sr^2; % sph scalar field
    laplacF = ee.getLaplacianSphS(F);
    %fprintf('\nnamb2 of F : %s\n', (simplify(laplacF)));
    pretty(simplify(laplacF)); % very pretty
    
    testVcyn = [ cr^3, (cos(cf))^6, cz*cr];  % Laplacians on vector fields are an algebra abortion
    test = ee.getLaplacianCynV(testVcyn);
    display(test);
    testVsph = [ sr*cos(sf), st*sr^4, st*sf*sr];
    test = ee.getLaplacianSphV(testVsph);
    display(test);
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    kS = 2;
    kV = [ kS, kS, kS ];
    rv = [ rx, ry, rz ];
    const = 3;
    
    Av = [ kS, 0, 0 ]; % vanishing divg and curl ,   divg(Av) = 0 , curl(Av) = [ 0, 0, 0 ]
    %fun_graphVF(Av, [0,0,2], 2);
    divgA = ee.getDivgRec(Av);
    curlA = ee.getCurlRec(Av);
    fprintf('vecFld A  IS solenoidal beause divergence is 0\n');
    fprintf('divg(A) = %.2f\n', divgA);
    fprintf('curl(A) = [ %.2f , %.2f, %.2f ]\n', curlA);
    
    Bv = kS .* rv;                % the radial vector is just componets to directions
    %fun_graphVF(Bv, [0,0,2], 2);
    divgB = ee.getDivgRec(Bv);
    curlB = ee.getCurlRec(Bv);
    fprintf('\nvecFld B  NOT solenoidal beause divergence is not 0\n');
    fprintf('divg(B) = %.2f\n', divgB);
    fprintf('curl(B) = [ %.2f , %.2f, %.2f ]\n', curlB);
    
    Cv = cross(kV, rv);
    %fun_graphVF(Cv, [0,0,2], 2);
    divgC = ee.getDivgRec(Cv);
    curlC = ee.getCurlRec(Cv);
    fprintf('\nvecFld C  IS solenoidal beause divergence is 0\n');
    fprintf('divg(C) = %.2f\n', divgC);
    fprintf('curl(C) = [ %.2f , %.2f, %.2f ]\n', curlC);
    
    Dv =  Cv + const .* rv ;
    %fun_graphVF(Dv, [0,0,2], 2);
    divgD = ee.getDivgRec(Dv);
    curlD = ee.getCurlRec(Dv);
    fprintf('\nvecFld D  NOT solenoidal beause divergence is not 0\n');
    fprintf('divg(D) = %.2f\n', divgD);
    fprintf('curl(D) = [ %.2f , %.2f, %.2f ]\n', curlD);  
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    Bv = [ry + rz * cos(rx*rz), rx, rx * cos(rx*rz)]; % given...
    divgB = ee.getDivgRec(Bv);
    curlB = ee.getCurlRec(Bv);
    fprintf(' given vector field B : [ %s, %s, %s ]\n', Bv);
    fprintf(' divg(B) = %s\n', divgB);
    fprintf('   since the divergence is not 0, the field is not solenoidal\n');
    fprintf(' curl(B) = [ %s, %s, %s ]\n', curlB);
    fprintf('   since the curl = 0, the field is irrational/potential/conservative\n');
    fun_graphVF(Bv, [0,0,2], 2);
end

%------------------------------------------------------------------------------------------ #23
if sel == 23
    
    
end
    
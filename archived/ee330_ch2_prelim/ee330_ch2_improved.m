%{
    let the globals ride, they will be helpful in the files
        let double() and simplify or double(simplify()) help you out

    #0  just some obj testing
    #1  ch2  pr2.6
    #2  ch2  pr2.7       
    #3  ch2  pr2.8        {nice graph}
    #4  ch2  p2.1
    #5  ch2  p2.2   just a point chek
    #6  ch2  p2.3   just a point check
    #7  ch2  p2.4   pts
    #8  ch2  p2.6   basic subs
    #9  ch2  p2.7    vec subs        the computer is your friend here
    #10 ch2  p2.8    vec subs
    #11 ch2  p2.9    vec subs  with point sub
    #12 ch2  p2.10   vec trans
    #13 ch2  p2.12   vec, pt, comp change
    #14 ch2  p2.13   vec trans
    #15 ch2  p2.14   unit vec relations
    #16 ch2  p2.16   CS SC point trans
    #17 ch2  p2.17    dp on units
    #18 ch2  p2.18    prove perp by cp
    #19 ch2  p2.19    cps/dps
    #20 ch2  p2.20    resolve comp from rec -> cyn
    #21 ch2  p2.21     transform rec vec field into cyn vector field
    #22 ch2  p2.22     let computer help
    #23 ch2  p2.23     ex2.2
    #24 ch2  p2.24     points and distance
    #25 ch2  p2.25    more pts...
    #26 ch2  p2.26    pt->conv
    #27 ch2  p2.27   more pts
    #28 ch2  p2.28   points...good one to see conv method
    #29 ch2  p2.29   points and comps
    #30 ch2  p2.30   a mixed vector field
    #31 ch2  p2.32    a vec field

%}
%color = uisetcolor([1, 1, 0], 'Selecf Color'); % .9, .9, .9 is nice
%clear all;  if things aren't going away
%consts = cls_CONST();  % common physical constants
%consts.check();
%consts.help();
clc;
close all;
clearvars;
con = cls_SUPERconvert();

sel = 31;

global rx; global ry; global rz; % rectangular params  rx, ry, rz        naugthy
global cr; global cf; global cz; % cylindrical params  cr, cf, cz
global sr; global st; global sf; % spherical params    sr, st, sf
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); syms cf; assume(cf,'real'); syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); syms st; assume(st,'real'); syms sf; assume(sf, 'real');


%------------------------------------------------------------------------------------------ #0
if sel == 0
    % the globals are working for a regular funtion or a static method in the class
    %{
    mtx = con.getUnits('SR');
    disp(mtx);
    answ = subs(mtx(2, 3), sf, pi/3);
    disp(answ);
    %con.feed(4, 0, pi/2, 'S');   % ?
    %con.feed(2^.5, 5*pi/4, 1, 'C');
    %con.feed(0, 3, 0, 'R');
    %con.feed(6, pi, 0, 'S');
    con.feed(4, 3*pi/2, -2, 'C');
    disp(con.inType);
    con.print();
    con.mag();
    con.graph();
    con.clear();
    con.print();
    %}
end   


%------------------------------------------------------------------------------------------ #1
if sel == 1
    Hv = [4, -3, 5]; % given a const vec H  in cyn
    Pt = [1, pi/2, 0]; % given a point P in cyn
    con.feed(Pt(1), Pt(2), Pt(3), 'C');
    con.graph();
    % the componet of H parallel to surface cr = 1   aka unit vector <1, 0, 0> in cyn
    %  ... anything not part of acr .... < 0, -3 * acf, 5 * acz >
    % verify with cp       cross product is the key  parallel = CP
    acr = [ 1, 0 , 0 ];
    answ = cross(Hv, acr);
    fprintf('answer for p2.6 is < %d, %d, %d >\n', answ);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    Gv = [ 20, 50, 50 ];    % given a const vector G in shp
    pt = [ 1, pi/2, pi/6];  % given a point in sph where vector is
    con.feed(pt(1), pt(2), pt(3), 'S');
    con.graph();
    % want comp G perp to surface st = pi/2  ....   perp = DP
    ast = [ 0, 1, 0];
    perpMag = dot(Gv, ast);
    perp = perpMag .* ast;
    fprintf('answer for p2.7 is < %d , %d, %d >\n', perp);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % if cr = 2 and cr{rz} = 1 ....  
    % a cylinder of radius 2 is extended on the z axis    [ cr eqn ]
    % a plane at z = 1 is paralle to the xy plane ..... so a cirle is cut from cylinder
    figure();
    maxP = 7;
    hold on;
    grid on;
    view(145, 30);
    title('triangle by 3 points as vectors', 'fontsize', 16);
    xlabel('x axis');
    ylabel('y axis');
    zlabel('z axis');
    xlim([-maxP, maxP]);
    ylim([-maxP, maxP]);
    zlim([-maxP, maxP]);
    xax = linspace(-maxP, maxP, 128);
    yax = linspace(-maxP, maxP, 128);
    zax = linspace(-maxP, maxP, 128);
    plot3(xax, 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax, 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 1);
    plot3(maxP,0,0,'yo', 'markersize', 6, 'linewidth', 3); % +x direction
    plot3(0,maxP,0,'yo', 'markersize', 6, 'linewidth', 3); % +y direction
    plot3(0,0,maxP,'yo', 'markersize', 6, 'linewidth', 3); % +z direction
    
    plot3(2, 0, 1, 'g.', 'markersize', 20);
    plot3(-2, 0, 1, 'g.', 'markersize', 20);
    plot3(0, 2, 1, 'g.', 'markersize', 20);
    plot3(0, -2, 1, 'g.', 'markersize', 20);
    
    cordX = linspace(-maxP, maxP, 100); 
    cordY = linspace(-maxP, maxP, 100);
    cordZ = linspace(-maxP, maxP, 100);
    x2s = zeros(100, 100) + 1;
    
    [x, y] = meshgrid(cordX, cordY);
    z1s = zeros(100, 100) + 1;
    surf(x, y, z1s, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none'); % z = 1
    
    r = 2;
    ht = maxP;
    [cordX, cordY, cordZ] = cylinder(r);
    cordZ = cordZ * ht;
    surf(cordX, cordY, cordZ,'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none');
    cordZ = cordZ*-1;
    surf(cordX, cordY, cordZ,'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none');
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    pp = [2, 5, 1]; % 3 rec pts  P, Q, R  to cyn and sph
    qp = [-3, 4, 0]; 
    rp = [6, 2, -4];
    
    con.feed(pp, 'R');
    con.print();
    
    con.feed(qp, 'R');
    con.print();
    
    con.feed(rp, 'R');
    con.print();
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    p1 = [2, pi/6, 5];
    con.feed(p1, 'C');
    con.print();
    p2 = [1, pi/2, -3];
    con.feed(p2, 'C');
    con.print();
    p3 = [10, pi/4, pi/3];
    con.feed(p3, 'S');
    con.print;
    p4 = [4, pi/6, pi/3];
    con.feed(p4, 'S');
    con.print();
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    pt = [2, 6, -4];
    con.feed(pt, 'R');
    con.print();
end
 

%------------------------------------------------------------------------------------------ #7
if sel == 7
    pt = [5, deg2rad(120), 1];
    con.feed(pt, 'C');
    con.print();
    pt = [10, pi/3, pi/6];
    con.feed(pt, 'S');
    con.print();
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    Vv = rx * rz - rx * ry + ry * rz;           % simple, but effective
    frx = cr * cos(cf);
    fry = cr * sin(cf);
    frz = cz;
    Vvsph = subs(Vv, [rx, ry, rz], [frx, fry, frz]);
    fprintf('in cyn: %s \n', Vvsph);
    % or use obj
    answ = con.simpleRCsub(Vv);
    fprintf('in cyn: %s \n', answ);
    
    Uv = rx^2 + 2*ry^2 + 3*rz^2;
    frx = sr * sin(st) * cos(sf);
    fry = sr * sin(st) * sin(sf);
    frz = sr * cos(st);
    Uvsph = subs(Uv, [rx, ry, rz], [frx, fry, frz]);
    fprintf('in sph: %s \n', Uvsph);
    % or use obj
    answ = con.simpleRSsub(Uv);
    fprintf('in sph: %s \n', answ);
end
  

%------------------------------------------------------------------------------------------ #9
if sel == 9
    Fv = (1 / sqrt(rx^2 + ry^2 + rz^2) ) * [rx, ry, 4];
    answ = con.vecRCsub(Fv);
    fprintf('in cyn: %s \n', simplify(answ));
    answ = con.vecRSsub(Fv);
    fprintf('\n');
    fprintf('in sph: %s \n', simplify(answ));
    
    fprintf('\n\n-------------------------------------------\n');
    
    Gv = ( (rx^2 + ry^2) / sqrt(rx^2 + ry^2 + rz^2) ) * [rx, ry, rz];
    answ = con.vecRCsub(Gv);
    fprintf('in cyn: %s \n', simplify(answ));
    answ = con.vecRSsub(Gv);
    fprintf('\n');
    fprintf('in sph: %s \n', simplify(answ));
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    Bv = [ sqrt(rx^2 + ry^2) , ry / sqrt(rx^2 + ry^2) , rz ];
    answ = con.vecRCsub(Bv);
    fprintf('in cyn: %s\n', simplify(answ));  % I think that is it...no answer 
    answ = con.vecRSsub(Bv);
    fprintf('in sph: %s\n', simplify(answ));
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    Av = [2, 3, 4];       % const cyn vec
    pt = [2, pi/2, -1];   % it wants this cyn point evaluated, but in rec
    
    mtx = con.getU_CR();
    sol = mtx * transpose(Av);
    answ = subs(sol, [cr, cf, cz], pt);
    display(simplify(answ));
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    Av = [cr * sin(cf), cr * cos(cf), -2*cz]; % cyn -> rec
    trans = con.getU_CR();
    temp = trans * transpose(Av);
    answ = simplify(temp);
    display(answ);
    trann = con.compCR(answ);
    display(trann);
    
    Bv = [4*sr*cos(sf), sr, 0];
    trans = con.getU_SR();
    
    temp = trans * transpose(Bv);
    answ = simplify(temp);
    display(answ);
    trann = con.compSR(answ);
    display(simplify(trann));
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    pt = [2, pi/2, 3*pi/2];  % sph pt given
    Bv = [sr*sin(st), 0, -sr^2 * cos(sf)];
    Bvp = subs(Bv, [sr, st, sf], pt);
    display(Bvp);
    con.feed(Bvp, 'S');
    con.print();
end

%------------------------------------------------------------------------------------------ #14
if sel == 14
    Bv = [ 0, 0, rx]; % given in rec
    tran = con.getU_RC;
    cyn = tran * transpose(Bv);
    display(cyn);
    cyn = con.compRC(cyn);
    fprintf('cyn: %s\n', cyn);
end

%------------------------------------------------------------------------------------------ #15
if sel == 15
    ar = con.getU_RC();
    arx = ar(1,:);
    acr = [1, 0, 0];
    arxDPacr = dot(arx, acr);
    display(arxDPacr);
    
    acf = [0, 1, 0];
    arxDPacf = dot(arx, acf);
    display(arxDPacf);
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    testP = [1, 2, -3];
    con.feed(testP, 'R');
    testPc = con.pts(2,:);
    testPs = con.pts(3,:);
    con.print();
    
    cyn = [cr, cf, cz];
    tempCS = con.simpleCSsub(cyn);
    ptCS = subs(tempCS, [sr, st, sf], testPs);
    fprintf('\nCS: ( %.4f, %.4f°, %.4f ) ...match\n',...
        ptCS(1), rad2deg(ptCS(2)), ptCS(3) );
    
    sph = [sr, st, sf];
    tempSC = con.simpleSCsub(sph);
    ptSC = subs(tempSC, [cr, cf, cz], testPc);
    fprintf('\nSC: ( %.4f, %.4f°, %.4f° ) ...match\n',...
        ptSC(1), rad2deg(ptSC(2)), rad2deg(ptSC(3)) );
    
    CSu = con.getU_CS();
    Ac = [1, 1, 1];
    display(CSu);
    tran = CSu * transpose(Ac);
    display(tran);
    subSC = con.simpleCSsub(tran);
    display(subSC);
    temp = subs(subSC, [cr, cf, cz], Ac);
    display(temp);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    acr = [1, 0, 0];
    arx = [1, 0, 0];  % want dp, but need to convert unit vec first
    pt = [2, 0, -1];
    con.feed(pt, 'R');
    tran = con.getU_RC();
    display(tran); % rec in terms of cyn
    arxT = tran * transpose(arx);
    arxT = transpose(arxT);
    display(arxT); % arx is now in terms of cyn
    result = dot(acr, arxT);
    display(result); % now just use the point or fact cos(cf) = rx / sqrt(rx^2 + ry^2)
    % ...it = 1, they appear to be pointing in the same direction at this particular point...think about it
    con.graph();
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    Av = [cr*sin(cf), cr*cos(cf), cr]; % a given cyn vec A
    Bv = [cr*sin(cf), cr*cos(cf), -cr]; % a given cyn vec B
    result = dot(Av, Bv);
    display(result);
    % take out the cr^2, 1 - 1 = 0, they are perp at all pts
    display(simplify(result));
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    Av = [3, 2, 1]; % a given const cyn vec A
    Bv = [5, 0, -8]; % a given const cyn vec B
    resA = Av + Bv;
    display(resA);
    resB = dot(Av, Bv);
    display(resB);
    resC = cross(Av, Bv);
    display(resC);
    magA = norm(Av);
    magB = norm(Bv);
    quo = resB / ( magA * magB);
    resD = acosd(quo);
    display(resD);
    quo = norm(resC) / (magA * magB);
    resD = asind(quo);
    display(resD);
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    Gv = [3*cr, cr*cos(cf), -cz^2]; % a given vector field G in cyn
    ptQr = [3, -4, 6]; % given point in rec
    con.feed(ptQr, 'R');
    ptQc = con.pts(2,:);
    Gvec = subs(Gv, [cr, cf, cz], ptQc); % Gv evaluated at point, it's cyn representation
    display(double(simplify(Gvec)));
    tran = con.getU_CR();
    display(tran);
    answ = tran * transpose(Gvec);
    display(answ); % arx is top row
    temp = answ(1,:);
    final = subs(temp,[cr, cf, cz], ptQc);
    display(double(simplify(final))); % component along arx = [1, 0, 0] at pt
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    Gv = [ry*rz, rx*rz, rx*ry]; % given rec vec G, trans to cyn
    tran = con.getU_RC();
    display(tran);
    res = tran * transpose(Gv);
    display(res);
    result = con.simpleRCsub(res);
    display(result);
    display(simplify(result));
end

%------------------------------------------------------------------------------------------ #22
if sel == 22
    % complete cyn to rec trans
    tran = con.getU_CR();
    syms Acr; syms Arx;
    syms Acf; syms Ary;
    syms Acz; syms Arz;
    Ac = [Acr, Acf, Acz];
    Ar = [Arx, Ary, Arz];
    map = tran * transpose(Ac);
    mapO = con.simpleCRsub(map);
    fprintf('Arx = %s \n', mapO(1,1) );
    fprintf('Ary = %s \n', mapO(2,1) );
    fprintf('Arz = %s \n', mapO(3,1) );
    % send it back just to check
    tran = con.getU_RC();
    map = tran * mapO;
    mapZ = con.simpleRCsub(map);
    mapZ = simplify(mapZ);
    fprintf('\nAcr = %s \n', mapZ(1,1) );
    fprintf('Acf = %s \n', mapZ(2,1) );
    fprintf('Acz = %s \n', mapZ(3,1) );  % it made it full circle
    
    % try with a sph -> rec -> sph   this time with a const vector instead of a symbol
    %{
    pt = [ 3, pi/2, pi/2];
    tran = con.getU_SR();
    trans = tran * transpose(pt);
    con.feed(pt, 'S');
    con.print();
    con.graph();
   
    check = subs(trans, [sr, st, sf], pt);
    check = double(check);
    display(check);
    %}
    pt = [2, 0, 0];
    con.feed(pt, 'R');
    con.print();
    tran = con.getU_RS();
    trans = tran * transpose(pt);
    tansx = con.simpleSRsub(trans);
    display(tansx);
    display(double(subs(tansx,[rx, ry, rz],pt)));  
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    
    ptA = [10, pi/2, 3*pi/4]; % sph
    Av = [cr*cz*sin(cf), 3*cr*cos(cf), cr*cos(cf)*sin(cf)]; % Av in cyn -> sph
    trans = con.getU_CS();
    %display(trans);
    tranA = trans * transpose(Av);
    %display(tranA);
    tranA = con.simpleCSsub(tranA);
    %display(tranA);
    tranA = subs(tranA, [sr, st, sf], ptA);
    display(double(tranA));
    
    ptB = [2, pi/6, 1]; % cyn
    Bv = [sr^2, 0, sin(st)]; % sph -> cyn
    trans = con.getU_SC();
    %display(trans);
    tranB = trans * transpose(Bv);
    %display(tranB);
    tranB = con.simpleSCsub(tranB);
    %display(simplify(tranB));
    tranB = subs(tranB, [cr, cf, cz], ptB);
    display(double(tranB));
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    a1 = [2, 1, 5];
    a2 = [6, -1, 2];
    a = sqrt(sum((a2 - a1).^2));
    b1 = [3, pi/2, -1];
    con.feed(b1, 'C');
    temp1 = con.pts(1,:); 
    b2 = [5, 3*pi/2, 5];
    con.feed(b2, 'C');
    temp2 = con.pts(1,:);
    b = sqrt(sum((temp1 - temp2).^2));
    display(a);
    display(b);
    
    p1 = [10, pi/4, 3*pi/4];
    p2 = [5, pi/6, 7*pi/4];
    con.feed(p1,'S');
    temp1 = con.pts(1,:);
    con.feed(p2, 'S');
    temp2 = con.pts(1,:);
    ptd = temp1-temp2;
    ptd = ptd.^2;
    ptd = sum(ptd);
    ptd = sqrt(ptd);
    display(ptd);
end


%------------------------------------------------------------------------------------------ #25
if sel == 25
    pt = [4, pi/6, 0];  % book is wrong
    qt = [6, pi/2, pi];
    con.feed(pt,'S');
    con.graph();
    temP = con.pts(1,:);
    con.feed(qt,'S');
    con.graph();
    temQ = con.pts(1,:);
    dist = temP - temQ;
    dist = dist.^2;
    dist = sum(dist);
    dist = sqrt(dist);
    display(dist);
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    acr = [1, 0, 0];
    acf = [0, 1, 0];
    inp = [0, 4, -1];
    tran = con.getU_CR();
    display(tran);
    mtx_acr = tran * transpose(acr);
    mtx_acf = tran * transpose(acf);
    display(mtx_acr);
    display(mtx_acf);
    mtx_acr = con.simpleCRsub(mtx_acr);
    mtx_acf = con.simpleCRsub(mtx_acf);
    display(mtx_acr);
    display(mtx_acf);
    ans_acr = subs(mtx_acr, [rx, ry, rz], inp);
    ans_acf = subs(mtx_acf, [rx, ry, rz], inp);
    display(ans_acr);
    display(ans_acf);
end


%------------------------------------------------------------------------------------------ #27
if sel == 27
    Av = [2*cz - sin(cf), 4*cr + 2*cos(cf), -3*cr*cz];
    Bv = [cr*cos(cf), sin(cf), 1];
    pt = [1, pi/3, -1];
    
    Avs = double(subs(Av, [cr, cf, cz], pt));
    Bvs = double(subs(Bv, [cr, cf, cz], pt));
    
    AdpB = dot(Avs, Bvs);
    display(AdpB);
    Amag = norm(Avs);
    Bmag = norm(Bvs);
    %display(Bmag);
    quo = AdpB/(Amag * Bmag);
    angl = acosd(quo);
    display(angl);
    AcpB = cross(Avs, Bvs);
    display(AcpB);
    
    feeler = cross(Av, Bv);
    feeler = double(subs(feeler, [cr, cf, cz], pt));
    display(feeler);
    display(feeler/norm(feeler));
    
    %{
    AdpB = dot(Av, Bv);
    AcpB = cross(Av, Bv); % technically you have the negative of it also....still normal
    display(simplify(AcpB));
    
    AdpBs = subs(AdpB, [cr, cf, cz], pt);
    AcpBs = subs(AcpB, [cr, cf, cz], pt);
    display(double(AcpBs));
    Amag = norm(subs(Av, [cr, cf, cz], pt));
    Bmag = norm(subs(Bv, [cr, cf, cz], pt));
    
    quo = AdpBs / (Amag * Bmag);
    ang = acosd(quo);
    display(double(ang));
    
    normu = norm(AcpBs);
    temp = AcpBs/normu;
    display(double(temp));
    
    
    Avs = subs(Av, [cr, cf, cz], pt);
    Amag = norm(Avs);
    Bvs = subs(Bv, [cr, cf, cz], pt);
    Bmag = norm(Bvs);
    
    quo = AdpB / (Amag * Bmag);
    angl = acosd(quo);
    display(double(simplify(angl)));
    
    quo = norm(AcpB) / (Amag * Bmag);
    angl = asind(quo);
    display(double(simplify(angl)));
    
    unorm = AcpB ./ norm(AcpB);
    display(double(simplify(unorm)));
    %}
    
end
   
%------------------------------------------------------------------------------------------ #28
if sel == 28
    A = [2, 4, 10];  % rec
    B = [-5, 1, -3]; % cyn
    ptr = [0, 2, -5]; % rec
    con.feed(ptr, 'R');
    con.print();
    %con.graph();
    %ptc = con.pts(2,:);
    % take B in cyn -> rec
    
    tran = con.getU_CR();
    display(tran);
    map = tran * transpose(B);
    map = con.simpleCRsub(map);
    map = subs(map, [rx, ry, rz], ptr);
    map = transpose(map);
    display(map);
    
    vec = A + map;
    display(vec);
    
    Bmag = norm(map);
    Amag = norm(A);
    AdpB = dot(A, map);
    quo = double(AdpB / (Amag * Bmag));
    angl = double(acosd(quo));
    display(180-angl); % it is in q1 ?   bullshit
    % A on B scalar comp
    unitBr = map ./ norm(map);
    comp = double(dot(A, unitBr));
    display(comp);
end


%------------------------------------------------------------------------------------------ #29
if sel == 29
    Bv = [cr^2 * sin(cf), (cz - 1) * cos(cf), cr^2]; % given vector in cyn
    pt = [ 4, pi/4, -1]; % given cyn pt
    % need to find dot(B, arx)    ... cyn -> rec, dp then
    Bvr = con.getU_CR() * transpose(Bv);
    Bvr = subs(Bvr, [cr, cf, cz], pt);
    Bvr = transpose(Bvr);
    display(Bvr);
    % not just a DP
    display( dot(Bvr, [1, 0, 0])); % 9 is good
end


%------------------------------------------------------------------------------------------ #30
if sel == 30
    Gv = [rx * cos(sf) / cr, 2*ry*rz/ cr^2, 1-(rx^2 / rz^2)]; % in rec, want sph
    
    Gv = con.simpleRSsub(Gv); % sub out rec parts
    Gv = con.simpleCSsub(Gv); % sub out cyn parts
    
    Gvs = con.getU_RS() * transpose(Gv); % it went to sph
    
    %Gvs = con.simpleRSsub(Gvs); % sub out rec parts
    %Gvs = con.simpleCSsub(Gvs); % sub out cyn parts
    
    Gvs = simplify(Gvs);
    display(Gvs);
end


%------------------------------------------------------------------------------------------ #31
if sel == 31
    Jv = [sr*sin(st)*cos(sf), -cos(2*st)*sin(sf), tan(st/2)*log(sr)]; % given spherical vector
    pts = [2, pi/2, 3*pi/2]; % given sph pt
    con.feed(pts, 'S');
    %con.graph();
    
    Jvps = double(subs(Jv, [sr, st, sf], pts)); % evaluate in sph
    %display(Jvps);
    Jvt = con.getU_SR() * transpose(Jvps); % putting it in rec
    Jvt = subs(Jvt, [sr, st, sf], pts); % values in
    
    % if J is paralel to raz, then it can only be the z comp
    fprintf('the part of J parallel to raz is: %.4f \n', Jvt(3));
    map = con.getU_RS() * Jvt;
    %display(map);
    
end
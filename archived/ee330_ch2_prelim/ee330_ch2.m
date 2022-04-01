global tester;
syms tester;
assume(tester, 'real');

global rx;
global ry;
global rz;
global cr;
global cf;
global cz;
global sr;
global st;
global sf;

syms rx;                % rectangular params  rx, ry, rz
assume(rx, 'real');
syms ry;
assume(ry,'real');
syms rz;
assume(rz, 'real');
syms cr;                % cylindrical params  cr, cf, cz
assume(cr, 'real');
syms cf;
assume(cf,'real');
syms cz;
assume(cz, 'real');
syms sr;                % spherical params  sr, st, sf
assume(sr, 'real');
syms st;
assume(st,'real');
syms sf;
assume(sf, 'real');
    
    
%{
    atan2()  ... the 4 quadrant solution

    argnames(fun)  gives variables to function

    #1  ch2  ex  1     rec points и vector   to cyn и sph 
    #2  ch2  pp  1     rec points и vector   to cyn и sph 
    #3  ch2  ex  2     sph vec 
    #4  ch2  pp  2     cyn vec to rec
    #5  ch2  ex  3     working with vectors of same coord sys   book looks wrong  {plane graphs}
    #6  ch2  pp  3     vector field is given in cyn
    #7  ch2  ex  4     more resolving and trans
    #8  ch2  pp  4     resolving w/ tangents
            
%}

%color = uisetcolor([1, 1, 0], 'Selecf Color'); % .9, .9, .9 is nice
%clear all;  if things aren't going away
%consts = cls_CONST();  % common physical constants
%consts.check();
%consts.help();

clc;
close all;
clearvars;


ptCon = cls_PTSconv("give me points");  % R, C, S, feed(), ptCon.pts(i,j) for comp

%vecCon = cls_VECconv();    shit can

sel = 8;    % CHANGE to see desired calculation , see header for options
prompt = '\n any key to proceed :';
%str = input(prompt, 's');
%pause(1);
%clc;

syms rx;                % rectangular params  rx, ry, rz
assume(rx, 'real');
syms ry;
assume(ry,'real');
syms rz;
assume(rz, 'real');
syms cr;                % cylindrical params  cr, cf, cz
assume(cr, 'real');
syms cf;
assume(cf,'real');
syms cz;
assume(cz, 'real');
syms sr;                % spherical params  sr, st, sf
assume(sr, 'real');
syms st;
assume(st,'real');
syms sf;
assume(sf, 'real');


%----------------------------------------------------------------------------------------------- #1
if sel == 1
    
    Prx = -2;   % the given point      starting from a given rectangular point P and vector A
    Pry = 6;
    Prz = 3; 
    ptCon.feed(Prx, Pry, Prz, 'R');
    
    Arx = ry;        % given x componenet of vector
    Ary = (rx + rz);  % given y component of vector
    Arz = 0;        % given z component of vector
    
    [Pcf, Pcr, Pcz] = cart2pol(Prx, Pry, Prz);                        
    fprintf('\nrect to cyl  for point\n');
    fprintf('MATLAB:  rect ( %d, %d, %d ), cylin ( %d, %d, %d )\n',...   % CORRECT
        Prx, Pry, Prz, Pcr, rad2deg(Pcf), Pcz);
    Pcr = sqrt(Prx^2 + Pry^2);
    Pcf = atan(Pry / Prx);
    Pcz = Prz;
    fprintf('Book  :  rect ( %d, %d, %d ), cylin ( %d, %d, %d )\n',...% must consider quad2
        Prx, Pry, Prz, Pcr, rad2deg(Pcf), Pcz);                       
    fprintf('answer (6.32 , 108.43° , 3)\n');
    
    [Psf, Pst, Psr] = cart2sph(Prx, Pry, Prz);
    fprintf('\nrect to sph  for point\n');
    fprintf('MATLAB:  rect ( %d, %d, %d ), sph ( %d, %d, %d )\n',...
        Prx, Pry, Prz, Psr, rad2deg(Pst), rad2deg(Psf) );
    Psf = atan(Pry / Prx);
    Pst = atan( sqrt( Prx^2 + Pry^2 ) / Prz);
    Psr = sqrt(Prx^2 + Pry^2 + Prz^2);
    fprintf('BOOK:  rect ( %d, %d, %d ), sph ( %d, %d, %d )\n',...
        Prx, Pry, Prz, Psr, rad2deg(Pst), rad2deg(Psf) );
    fprintf('answer ( 7, 64.62°, 108.43° )\n');   
    
    ptCon.print();
    feed = [ Arx, Ary, Arz ];
    rV = fun_vecCon(ptCon.pts, feed, 'R');
    display(rV);
    fprintf('\nlength of all 3 vectors should be same :\n');
    fprintf('length of rectangular: %.4f\n', sqrt( rV(1, 1)^2 + rV(1, 2)^2 + rV(1, 3)^2 ) );
    fprintf('length of cylindrical: %.4f\n', sqrt( rV(2, 1)^2 + rV(2, 2)^2 + rV(2, 3)^2 ) );
    fprintf('length of cylindrical: %.4f\n', sqrt( rV(3, 1)^2 + rV(3, 2)^2 + rV(3, 3)^2 ) );
end


%----------------------------------------------------------------------------------------------- #2
if sel == 2
    
    ptCon.clear();
    prx = 1; pry = 3; prz = 5;
    ptCon.feed(prx, pry, prz, 'R');
    fprintf('\npoint P is: \n');
    ptCon.print();
    %ptCon.magD();   the magnitudes are not at all the same...wtf ?
    
    str = input(prompt, 's');
    pause(1);
    clc;
    
    ptCon.clear();
    Trx = 0; Try = -4; Trz = 3;
    ptCon.feed(Trx, Try, Trz, 'R');
    fprintf('\npoint T is: \n');
    ptCon.print();
    
    str = input(prompt, 's');
    pause(1);
    clc;
    
    ptCon.clear();
    srx = -3; sry = -4; srz = -10;
    ptCon.feed(srx, sry, srz, 'R');
    fprintf('\npoint S is: \n');
    ptCon.print();
    
    str = input(prompt, 's');
    pause(1);
    clc;
    
    ptCon.clear();
    Trx = 0; Try = -4; Trz = 3;
    ptCon.feed(Trx, Try, Trz, 'R');
    fprintf('\npoint T is: \n');
    ptCon.print();
    
    Qrx = ( sqrt( rx^2 + ry^2 ) / sqrt ( rx^2 + ry^2 + rz^2 ) );
    Qry = 0;
    Qrz = - ( ry * rz ) / sqrt ( rx^2 + ry^2 + rz^2 ) ;
    feed = [ Qrx, Qry, Qrz ];
    vecResult = fun_vecCon(ptCon.pts, feed, 'R');
end


%----------------------------------------------------------------------------------------------- #3
if sel == 3
    ptCon.clear();
    cd1 = -3; cd2 = 4; cd3 = 0;
    ptCon.feed(cd1, cd2, cd3, 'R');
    bsr = 10 / sr; bst = sr * cos(st); bsf = 1;
    vec = [ bsr, bst, bsf ];
    ptCon.print();
    vecR = fun_vecCon(ptCon.pts, vec, 'S');
    
    resR = subs(bsr, [sr, st, sf], [ptCon.pts(3, 1), ptCon.pts(3, 2), ptCon.pts(3, 3)] );
    resE = subs(bst, [sr, st, sf], [ptCon.pts(3, 1), ptCon.pts(3, 2), ptCon.pts(3, 3)] );
    resA = subs(bsf, [sr, st, sf], [ptCon.pts(3, 1), ptCon.pts(3, 2), ptCon.pts(3, 3)] );
    mag = sqrt(resR^2 + resE^2 + resA^2);
    fprintf('mag should be %.4f \n', mag); % book is wrong
end

%----------------------------------------------------------------------------------------------- #4
if sel == 4
    Acr = cr * cz * sin(cf); Acf = 3 * cr * cos(cf); Acz = cr * cos(cf) * sin(cf);
    vecin = [ Acr, Acf, Acz ];
    % making this cyn vector  a rec vector ... get simplified output in comps, trans, part, и jac
    [ outVector, trans33, partial33, jak ] = this2that(vecin, 'CR'); 
    display(outVector);
    %display(trans33);
    %display(partial33);
    %display(jak);
    ptCon.feed(-3, 4, 0, 'R');
    ptCon.print();
    
    pause(1);
    clc;
    
    Asr = sr^2; Ast = 0; Asf = sin(st);
    vecin = [ Asr, Ast, Asf ];
    % making this cyn vector  a rec vector ... get simplified output in comps, trans, part, и jac
    [ outVector, trans33, partial33, jak ] = this2that(vecin, 'SR'); 
    display(outVector);
    display(trans33);
    display(partial33);
    display(jak);
end

%----------------------------------------------------------------------------------------------- #5
if sel == 5
    % given two uniform vectors in cynlindrical
    Ev = [-5, 10, 3];
    Fv = [1, 2, -6];
    EFmag = norm( cross(Ev, Fv) );
    fprintf('\nthe magnitude of E x F is %.4f\n', EFmag);
    
    % the component of E at (5, pi/2, 3) parallel to x=2, z=3
    ptCon.clear();
    ptCon.feed(5, pi/2, 3, 'C');
    ptCon.print(); % smart feature, stops all this....
    %fprintf('the cyn point ( %.4f, %.4f, %.4f ) \n', ptCon.pts(2,:) );
    %fprintf('.... is the rec point ( %.4f, %.4f, %.4f ) \n', ptCon.pts(1,:) );
    
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
    plot3(maxP-1,0,0,'yo', 'markersize', 6, 'linewidth', 3); % +x direction
    plot3(0,maxP-1,0,'yo', 'markersize', 6, 'linewidth', 3); % +y direction
    plot3(0,0,maxP-1,'yo', 'markersize', 6, 'linewidth', 3); % +z direction
    
    plot3(ptCon.pts(1, 1), ptCon.pts(1, 2), ptCon.pts(1, 3), 'g.', 'markersize', 20);
    plot3(2, 0, 0, 'k.', 'markersize', 20);
    plot3(0, 0, 3, 'k.', 'markersize', 20);
    plot3([0, ptCon.pts(1, 1)], [0, ptCon.pts(1, 2)], [0, ptCon.pts(1, 3)], 'g:', 'linewidth', 2);
    cordX = linspace(-maxP, maxP, 100); 
    cordY = linspace(-maxP, maxP, 100);
    cordZ = linspace(-maxP, maxP, 100);
    x2s = zeros(100, 100) + 2;
    [y, z] = meshgrid(cordY, cordZ);
    surf(x2s, y, z, 'FaceColor', 'b', 'FaceAlpha', .2, 'EdgeColor', 'none'); % x = 2
    [x, y] = meshgrid(cordX, cordY);
    z3s = zeros(100, 100) + 3;
    surf(x, y, z3s, 'FaceColor', 'r', 'FaceAlpha', .2, 'EdgeColor', 'none'); % x = 2
    
    % component of E paralles to line formed by plane intersections must be dot(E, ay) ay
    % but ay = < sin(cf) , cos(cf) , 0 >      comp is dot(Ev, ay) * ay
    % comp = dot(Ev, ay) ay  , but ay = < sin(cf) acr, cos(cf) act , 0 >
    % plug in point, = -5*1 acr  or -5 ary
    % since z axis is normal to surface z = 3  , dot it with comp E  and divide out for angle
    zHat = [0, 0, 1];
    Evec = ptCon.pts(1, :);
    E = norm(Evec);
    quo = dot(Evec, zHat) / ( E * 1 );
    ang = acosd(quo);
    disp(ang);
end

%----------------------------------------------------------------------------------------------- #6
if sel == 6
    % H given as a cyn vector
    Hcr = cr * cz * cos(cf);
    Hcf = exp(-2) * sin(cf/2);
    Hcz = cr^2;
    Hvec = [ Hcr, Hcf, Hcz ];
    ptCon.clear();
    ptCon.feed(1, pi/3, 0, 'C');
    ptCon.print();
    
    % arx is unit vector of x comp, but cyn to rec, arx = < cos(cf) , -sin(cf) , 0 >
    arx = [ cos(cf), -sin(cf), 0 ];
    H_dot_arx = dot(Hvec, arx);
    resultA = subs(H_dot_arx, [cr, cf, cz], [ ptCon.pts(2, :) ] );
    fprintf('H dot unit x : %.4f\n', resultA);
    
    %alternatively with fun_getUnits(type)
    arv = fun_getUnits('CR'); % get unit vectors of cyn in terms of rec
    arx = arv(1, :); % take unit 
    H_dp_arx = dot(Hvec, arx);
    display(double(subs(H_dot_arx, [cr, cf, cz], [ ptCon.pts(2, :) ] )));
    
    % part b is H x ast
    SRmtx = fun_getUnits('SC');
    uMat = SRmtx(2,:);
    display(uMat); % ast in terms of unit cyn
    result = cross(Hvec, uMat);
    display(result);
    % sc = atan(cr / cz)    but rho is 1 and z is 0, so sc = 90°     very important
    result = subs(result, st, pi/2);
    display(result);
    result = subs(result, [cr, cf, cz], [ptCon.pts(2, :)] );
    fprintf('answer is ( %.4f, %s, %s )\n', result);
    
    % part c, the vector component normal to rho = 1    .... dot(H, acz)
    % say you subed that point in to make a vector subV
    subV = subs(Hvec, [cr, cf, cz], ptCon.pts(2,:) );
    %display(subV);
    % now just dot it with < 1, 0, 0 >
    res = dot(subV, [1, 0, 0]);
    display(res);
    
    % scalar comp H tangent to z = 0 plane... ob that subV vector,
    % resolve it   
end


%----------------------------------------------------------------------------------------------- #7
if sel == 7
    Dsr = sr * sin(sf);
    Dst = (-1/sr) * sin(st) * cos(sf);
    Dsf = sr^2;
    Dvec = [ Dsr, Dst, Dsf ]; % vector field D as given
    Dp = [10, deg2rad(150), deg2rad(330)]; % point as given in sph
    ptCon.clear();
    ptCon.feed(Dp(1), Dp(2), Dp(3), 'S');
    ptCon.print;
    Dvp = subs(Dvec, [sr, st, sf], Dp); % D evaluated at point
    display(double(Dvp));
    
    % component of this Dvp tangent to spere r = 10 at point   
    % at r = 10 (given point) , only Dsr part is considered normal
    % V = Vtan + Vnorm   ->  Vtan = V - Vnorm
    fprintf('normal comp = %.1f * %s \n', Dvp(1), Dsr);
    
    % perp to D and tengent to cone theta = 150°  
    % same as D x ast
    display(cross(Dvp, [0,1,0]));
    display(double(cross(Dvp, [0,1,0])/norm(cross(Dvp, [0,1,0]))));
end

%----------------------------------------------------------------------------------------------- #8
if sel == 8
    Asv = [3, 2, -6];
    Bsv = [4, 0, 3]'
    AdpB = dot(Asv, Bsv);
    AcpB = cross(Asv, Bsv);
    display(AdpB);
    display(norm(AcpB));
    pt = [1, pi/3, 5*pi/4 ];
    unit = fun_getUnits('SR');
    zcp = unit(3,:);
    display(zcp);
    vex = subs(zcp, sf, pt(3));
    vexx = dot(Asv, zcp);
    vexxx = vexx * zcp;
    ans = subs(vexxx, st, pt(2));
    display(double(ans));
end
    
    
    
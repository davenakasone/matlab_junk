%{
    ee330 ch1

        remember these savages use ax = i hat, ay = j hat , az = k hat
    
        #0  appx  linear sys and basic solving
        #1  ch1, ex 1
        #2  ch1, pp 1
        #3  ch1, ex 2
        #4  ch1, pp 2
        #5  ch1, ex 3
        #6  ch1, pp 3
        #7  ch1, ex 4
        #8  ch1, pp 4
        #9  ch1, ex 5
        #10 ch1, pp 5
        #11 ch1, pp 6        graph protocol is getting there
        #12 ch1, ex 7
        #13 ch1, pp 7


    be careful with rad и  deg ....  sind() vs sin()
    cart2pol()  and pol2cart()  are going to be useful      automatic phasor info
    norm() is good for not having to do    length = sqrt(sum(V.^2))
    dot(A,B)   easier dot product than sum(A.*B)
    cross(A, B)  good if two 3-element vectors    ... no 7D allowed   can even do an entire system
    [x] = solve(x^2 + 2*x - 3 = 0, x)   will be helpful
    solving linear systems of equations is fundamental....
        sol = dsolve(eqn, 'Implict', true);  to get F(y(t)) = g(t) form
            lots of shit, one of the best functions, good series, systems, and approx

    dont forget .^  .*  ./   for element-element operations   w/o dot, entire operation occurs on mtx
    also, .'  transpose ect can change if complex numbers are in there

        protect your 1j for sqrt(-1)
        protect your e for sci notation 
        protect your pi  for sym(pi)     also all the regular function names

    assume() is a nice tool to change inputs    ... integer, odd, even, ect

%}

clc; 
clf;
close all;
clearvars;
select = 13; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice


%----------------------------------------------------------------------------------------------- 0
if select == 0
    
    mass_earth = 6.35e21; % how to use scientific notation
    fund_charge = 6.19e-19; % sci notation on negative
    display(mass_earth);
    display(6.35*(10^21));
    display(fund_charge);
    display(6.19*(10^-19));
    
    syms x;
    eqn = x^2 + 2*x -3;
    sols = solve(eqn,x);    % also use dsolve if there is a function and some of unknown derivatives
    display(sols);
    
    syms y(t);         % y is a function of t
    syms a;            %  a could be anything
    eqn = diff(y, t, 1) == a*y;    %  differential equation:  first derivative y w/r t   y' = a * y
    sol = dsolve(eqn);   % can even put condtions in
    display(sol);
    cond = y(0) == 5;
    sol = dsolve(eqn,cond);
    display(sol);
    
    syms b;   % 2 DE w/ conditions...
    eqn = diff(y, t, 2) == (a^2)*y;
    dy = diff(y,t,1);
    cond = [y(0)==b, dy(0)==1];
    ySol(t) = dsolve(eqn, cond);
    pretty(ySol);
    
    %{
        25x - 5y - 20z = 50
        -5x + 10y -4z = 0
        -5x -4y + 9z = 0
    %}
    A = [ 25, -5, -20; -5, 10, -4; -5, -4, 9];  % coeff matrix  RHS
    B = [50; 0; 0]; % values / LHS matrix
    % [m x n ] * [ n x m ] >  [ 3 x 3 ] * [ 3 x 1 ] = [ 3 x 1 ]   "rows = rows"
    % also A needs to be inverted...most be square and not singular 
    C = A\B;
    %C = inv(A) * B;    obsolete
    %C = (A^-1)*B;       also obsolete
    display(C);     % x = C(1,1) = C(1) = 29.6 , y = C(2,1) = C(2) = 26 , z = C(3,1) = C(3) = 28
end
    
    
%----------------------------------------------------------------------------------------------- 1
if select == 1
    
    Av = [10, -4, 6];
    A = sqrt(sum(Av.^2));   % how to get magnitude
    A1 = norm(Av);          % or use built in function
    av = Av ./ A;            % how to get a normal vector          
    
    
    Bv = [2, 1, 0];       % no z comp implies 0
    B = sqrt(sum(Bv.^2));
    B1 = norm(Bv);
    bv = Bv ./ B;
    
    
    fprintf('\nthe component of A on y-axis is %d\n', Av(2) );    % comp A on y axis
    
    Cv = 3*Av - Bv;   
    C = sqrt(sum(Cv.^2));
    C1 = norm(Cv);
    fprintf('the magnitude of 3A - B is %d \n', C);
    
    Dv = Av + 2*Bv;
    D = norm(Dv);
    dv = Dv ./ D;
    display(dv);
    display(norm(dv)); % length of 1 as expected
end


%----------------------------------------------------------------------------------------------- 2
if select == 2
    
    Av = [1, 0, 3];
    Bv = [5, 2, -6];
    
    Cv = Av + Bv;
    C = norm(Cv);
    fprintf('\nmag of A + B = %d \n', C);
    
    Cv = 5*Av - Bv;
    disp(Cv);
    
    fprintf('component of A in direction of y-axis = %d\n', Av(2) );
    
    Cv = 3*Av + Bv;
    C = norm(Cv);
    c = Cv ./ C;
    fprintf('unit vector parallel to 3A + B =\n');
    disp(c);
    disp(-c); % have to consider both directions
end


%----------------------------------------------------------------------------------------------- 3
if select == 3
    P = [0, 2, 4];
    Q = [-3, 1, 5];
    RP = P; % position vector is origin to point
    RQ = Q;
    
    fprintf('the position vector of point P(0, 2, 4) is\n');
    disp(RP);
    
    distv = RQ - RP;
    fprintf('the distance vector P to Q is\n');
    disp(distv);
    dist = norm(distv);
    fprintf('the distance from P to Q is %d \n', dist);
    
    Cv = 10 * (distv ./ dist);
    fprintf('the vector parallel to PQ with a mag of 10 is \n');
    disp(Cv);
    disp(-Cv); % they don't make a distinction between parallel and anit parallel...wtf
end


%----------------------------------------------------------------------------------------------- 4
if select == 4
    Pv = [1, -3, 5];
    Qv = [2, 4, 6];
    Rv = [0, 3, 8];
    
    fprintf('the position vector of P is\n');
    disp(Pv);
    fprintf('the position vector of R is\n');
    disp(Rv);
    
    distQRv = Rv - Qv; % 
    fprintf('the distance vector of Q to R is\n');
    disp(distQRv);
    distQR = norm(distQRv);
    fprintf('the distance from Q to R is\n');
    disp(distQR);
end


%----------------------------------------------------------------------------------------------- 5
if select == 5
    % river flows east @ 10km/h   ... <1,0,0> direction
    % guy walks up deck at 2km/h
    
    Bx = 10 * cosd(45);
    By = -10 * sind(45);
    Bv = [Bx, By];   % velocity of boat in km/h
    Gx = (-2) * cosd(45);
    Gy = (-2) * sind(45);
    Gv = [Gx, Gy];  % velocity of guy relative to boat
    
    Tv = Bv + Gv;   % combined velocity, w/r to earth
    T = norm(Tv);
    
    [theta, rho] = cart2pol(Tv(1),Tv(2));    %  BIG $$$
    
    degs = (180/pi)*theta;
    fprintf('the combined velocity is\n');
    disp(Tv);
    fprintf('the speed is %d    and the angle is %d   degrees\n', rho, degs);
end


%----------------------------------------------------------------------------------------------- 6
if select == 6
    % plane has a grond speed of 350 km/h west...    350 * < 1 , 0 >
    % wind blows north west at 40 km/h  .... 40 * < -cos(45) , sin(45) >
    Pv = 350 * [cosd(180), sind(180)]; % plane's velocity
    Wv = 40 * [cosd(135), sind(135)];  % wind's velocity vector
    Cv = Pv + Wv; % combined velocity
    [theta, rho] = cart2pol(Cv(1), Cv(2));
    fprintf('true heading is %d  degrees, at a speed of %d\n', (180/pi)*theta, rho);
    degs = 180-(180/pi)*theta;
    fprintf('aka %d degrees north of west\n', degs);
end


%----------------------------------------------------------------------------------------------- 7
if select == 7                  % use and abuse  cross() and dot()
    Av = [3, 4, 1];
    A = norm(Av);
    Bv = [0, 2, -5];
    B = norm(Bv);
    % want angle beteween these 2 ... take advantage of CP or DP properties
    
    ABdot = dot(Av, Bv);
    quo = ABdot / (A*B);
    syms th;
    sol = solve(cosd(th) == quo, th);  % take positive
    fprintf('the angle between A and B is %d  degrees\n', sol);
    % or
    sol = acosd(quo);
    fprintf('the angle between A and B is %d  degrees\n', sol);   % by-hand better
    
    ABcross = cross(Av, Bv);
    mag = norm(ABcross);
    quo = mag / (A * B);
    sol = asind(quo);
    fprintf('the angle between A and B is %d  degrees\n', sol);
end


%----------------------------------------------------------------------------------------------- 8
if select == 8
    Av = [1, 0, 3];
    A = norm(Av);
    Bv = [5, 2, -6];
    B = norm(Bv);
    % will find angle by CP and DP
    
    ABcp = cross(Av, Bv);
    mag = norm(ABcp);
    quo1 = mag / (A * B);
    sol1 = asind(quo1);
    fprintf('the angle between A and B is %d  degrees\n', 180-sol1);     %  need 180 -
    
    ABdp = dot(Av, Bv);    % trust DP more
    quo2 = ABdp / (A * B);
    sol2 = acosd(quo2);
    fprintf('the angle between A and B is %d  degrees\n', sol2);
end

%----------------------------------------------------------------------------------------------- 9
if select == 9
    Pv = [2, 0, -1];
    P = norm(Pv);
    Qv = [2, -1, 2];
    Q = norm(Qv);
    Rv = [2, -3, 1];
    R = norm(Rv);
    
    a = cross((Pv + Qv),(Pv-Qv));
    display(a);
    
    b = dot(Qv,cross(Rv,Pv));
    display(b);
    
    c = dot(Pv,cross(Qv,Rv));
    display(c);
    
    QRcp = cross(Qv,Rv);
    QRmag = norm(QRcp);
    quo = QRmag / (Q * R);
    display(quo);
    
    e = cross(Pv, cross(Qv, Rv));
    display(e);
    
    f = QRcp / norm(QRcp);
    display(f);  % technically +/- since no anti parallel
    display(norm(f));
    
    PQdp = dot(Pv, Qv);
    quo = PQdp / (P*Q);
    PQ = P*quo;
    fprintf('the component is %d\n', PQ);
    qv = Qv ./ Q;
    PQ = dot(qv, Pv);
    fprintf('the component is %d\n', PQ);
    PQv = qv .* PQ;
    display(PQv);
    display(norm(PQv));
end


%----------------------------------------------------------------------------------------------- 10
if select == 10
    Ev = [0, 3, 4];
    E = norm(Ev);
    Fv = [4, -10, 5];
    F = norm(Fv);
    
    fv = Fv ./ F; % E gets projected on F , need a unit vector in direction of F
    EF = dot(fv, Ev);     % the projected component of E on F has this magnitude
    EFv = EF .* fv;       % the vector that represents E on F
    display(EFv);
    display(norm(EFv));
    display(EF);          % pay attention to sign change
    
    EFcp = cross(Ev, Fv);
    mag = norm(EFcp);
    norm_perp = EFcp ./ mag;
    display(norm_perp);
    display(-norm_perp);  % consider both
end


%----------------------------------------------------------------------------------------------- 11
if select == 11
    
    maxP = 9;              % fit to point size
    aP = [4, 0, -1];       % points as given     here, the points were actually supposed to be vectors
    bP = [1, 3, 4];
    cP = [-5, -3, -3];
    
    av = bP - aP;  % side a is vector aP to bP
    bv = cP - bP;  % side b is vector bP to cP
    cv = aP - cP;  % side c is vector cP to aP   
    display(av + bv + cv); % closed loop formed, check   a + b + c = < 0, 0, 0 >
    display(aP + bP + cP);
    
    area = .5*norm(cross(aP, bP));
    display(area);
    area = .5*norm(cross(bP, cP));
    display(area);
    area = .5*norm(cross(cP, aP));
    display(area);
    
    %  90°  is at a to b
    ab_dp = dot(aP, bP);
    quo = ab_dp / (norm(aP) * norm(bP));
    th = acosd(quo);
    display(th);
    
    figure('Position', [20, 20, 700, 700]);
    hold on;
    view(145,30);
    grid on;
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
    % plot the given points
    aPstr = "( " + string(aP(1)) + " , " + string(aP(2)) + " , " + string(aP(3)) + " )";
    text(aP(1)+.3, aP(2)+.3, aP(3)+.3, aPstr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(aP(1), aP(2), aP(3), 'r.', 'markersize', 20);  
    bPstr = "( " + string(bP(1)) + " , " + string(bP(2)) + " , " + string(bP(3)) + " )";
    text(bP(1)+.3, bP(2)+.3, bP(3)+.3, bPstr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(bP(1), bP(2), bP(3), 'r.', 'markersize', 20);
    cPstr = "( " + string(cP(1)) + " , " + string(cP(2)) + " , " + string(cP(3)) + " )";
    text(cP(1)+.3, cP(2)+.3,cP(3)+.3, cPstr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(cP(1), cP(2), cP(3), 'r.', 'markersize', 20);
    % plot vectors forming sides of triangle
    plot3([aP(1), bP(1)], [aP(2), bP(2)], [aP(3), bP(3)], 'g'); % aP to bP
    plot3([bP(1), cP(1)], [bP(2), cP(2)], [bP(3), cP(3)], 'g'); % bP to cP
    plot3([cP(1), aP(1)], [cP(2), aP(2)], [cP(3), aP(3)], 'g'); % cP to aP
    % position vectors
    plot3([0, aP(1)], [0, aP(2)], [0, aP(3)], 'c'); % org to aP
    plot3([0, bP(1)], [0, bP(2)], [0, bP(3)], 'c'); % org to bP
    plot3([0, cP(1)], [0, cP(2)], [0, cP(3)], 'c'); % org to cP
    hold off;
end


%----------------------------------------------------------------------------------------------- 12
if select == 12
    maxP = 9;              % fit to point size
    p1 = [5, 2, -4];       % given some points, show they are on a line, and shortest path to p4
    p2 = [1, 1, 2];
    p3 = [-3, 0, 8];
    p4 = [3, -1, 0];
    
    % pick a point and get vectors to it
    Rp1p2 = p2 - p1; % the position vectors are subtracted, giving a vector p1 to p2
    Rp1p3 = p3 - p1; % the position vestors are subtracted, giving a vector p1 to p3
    Rp1p4 = p4 - p1; % the position vestors are subtracted, giving a vector p1 to p4
    
    % check cp's  if it is 0, they are paralell or anti parallel  (on same line)
    Rp1p2_cp_Rp1p3 = cross(Rp1p2, Rp1p3);
    display(Rp1p2_cp_Rp1p3);   % confirmed   cp = 0
    
    % you can even make an equation out of the line   start at p1, stop at p3 , t[0,1]
    pts = 20;  % change here
    linePts = zeros(pts, 3);
    step = 1/pts;
    t = step/2;
    for i = 1:pts
        
        temp = (1-t).*p1 + t.*p3;
        for k = 1:3
            linePts(i,k) = temp(k);
        end
        if i == 1
            t = 0;
        end
        t = t + step;
    end
    % or make a general eqn for p1 to p3 line      or p1 to p2
    syms lam real;
    lineEqn12 = p1 + lam.*(p2 - p1);
    display(lineEqn12);
    lineEqn13 = p1 + lam.*(p3 - p1);
    display(lineEqn13);
    display(subs(lineEqn12, lam, 2)); % if lam = 2,  p3 can be found on the linep1p2
    sol = solve(lineEqn12 == p3, lam);
    display(sol);                       % or just use solve()
    
    % the shortest distance from the line with p1, p2, p3  to p4 has to be the CP
    % just use the cross product of p1p4 x p1p2    or     p1p4 x p1p3    
    % it is a vector of points online formed with points online to p4
    % keep same reference point so you can make a triangle
    numer = norm(cross(Rp1p2, Rp1p4)); % this is the perpendicular distance
    denom = norm (Rp1p2);
    dist = numer/denom;
    display(dist);
    % or find the angle and do   dist = |p1p4| * sin(th)    th angle of p1p4 and p1p2
    p1p4_dp_p1p3 = dot(Rp1p4, Rp1p2);
    quo = p1p4_dp_p1p3 / ( norm(Rp1p4) * norm(Rp1p2) );
    th = acosd(quo);
    dist = norm(Rp1p4) * sind(th);  % old fashioned
    display(dist);
    
    % do get actual point, solve off distance formula, angle, or anything else the keeps relation
    closest = zeros(3);
    lm = solve(norm(lineEqn13) == dist, lam);
    closest = subs(lineEqn13, lam, lm);
    fprintf('closest point: %d, %d, %d \n', closest(1), closest(2), closest(3));
    lm = solve(norm(lineEqn12) == dist, lam);
    closest = subs(lineEqn12, lam, lm);
    fprintf('closest point: %d, %d, %d \n', closest(1), closest(2), closest(3));
    %lots of ways to solve
    
    figure('Position', [20, 20, 700, 700]);
    hold on;
    view(145,30);
    grid on;
    title('3 points and a path', 'fontsize', 16);
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
    % plot the given points
    p1str = "p1 ( " + string(p1(1)) + " , " + string(p1(2)) + " , " + string(p1(3)) + " )";
    text(p1(1)+.3, p1(2)+.3, p1(3)+.3, p1str,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(p1(1), p1(2), p1(3), 'k.', 'markersize', 20);  
    p2str = "p2 ( " + string(p2(1)) + " , " + string(p2(2)) + " , " + string(p2(3)) + " )";
    text(p2(1)+.3, p2(2)+.3, p2(3)+.3, p2str,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(p2(1), p2(2), p2(3), 'k.', 'markersize', 20);
    p3str = "p3 ( " + string(p3(1)) + " , " + string(p3(2)) + " , " + string(p3(3)) + " )";
    text(p3(1)+.3, p3(2)+.3,p3(3)+.3, p3str,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(p3(1), p3(2), p3(3), 'k.', 'markersize', 20);
    p4str = "p4 ( " + string(p4(1)) + " , " + string(p4(2)) + " , " + string(p4(3)) + " )";
    text(p4(1)+.3, p4(2)+.3,p4(3)+.3, p4str,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(p4(1), p4(2), p4(3), 'ro', 'markersize', 15, 'linewidth',2);
    % position vectors differences
    plot3([p1(1), p1(1) + Rp1p2(1)], [p1(2), p1(2) + Rp1p2(2)], [p1(3), p1(3) + Rp1p2(3)], 'c:', 'linewidth', 3); % p1 to p2
    plot3([p1(1), p1(1) + Rp1p3(1)], [p1(2), p1(2) + Rp1p3(2)], [p1(3), p1(3) + Rp1p3(3)], 'm--', 'linewidth', 1); % p1 to p2
    plot3([p1(1), p1(1) + Rp1p4(1)], [p1(2), p1(2) + Rp1p4(2)], [p1(3), p1(3) + Rp1p4(3)], 'g', 'linewidth', 1); % p1 to p4
    % plot parameterized line
    for i = 1:pts
        plot3([0, linePts(i,1)], [0, linePts(i,2)], [0, linePts(i,3)], 'k--', 'linewidth',2);
    end
    % plot closest path line to p4
    plot3([p4(1), closest(1)], [p4(2), closest(2)], [p4(3), closest(3)], 'r:', 'linewidth',4); % all points should be on this
    hold off;
end


%----------------------------------------------------------------------------------------------- 13
if select == 13
    maxP = 8;
    pa = [1, 2, -3];
    pb = [-4, 0, 5];
    pc = [7, -1, 2];
    pab = pb - pa; % vector pa to pb
    pac = pc - pa; % vector pa to pc
    
    fprintf('the distance of pa to pb is %d\n\n', norm(pab));
    
    syms lam;
    lineEqn = pa + lam .* (pb - pa);
    display('the line equation: ');
    pretty(lineEqn);
    
    %remember the various ways to get shortest distance
    pab_dp_pac = dot(pab, pac);
    quo = pab_dp_pac / ( norm(pab) * norm(pac) );
    th = acos(quo);
    dist = norm(pac) * sin(th);
    fprintf('the shortest distance is %d\n', dist);
    lm = solve(   sqrt(sum((lineEqn - pc).^2)) == dist, lam);
    close = double(subs(lineEqn, lam, lm));
    display(close);
    
    figure('Position', [800, 800, 700, 700]);
    hold on;
    view(145,30);
    grid on;
    title('3 points and a path', 'fontsize', 16);
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
    plot3(maxP,0,0,'y.', 'markersize', 50, 'linewidth', 13); % +x direction
    plot3(0,maxP,0,'y.', 'markersize', 50, 'linewidth', 13); % +y direction
    plot3(0,0,maxP,'y.', 'markersize', 50, 'linewidth', 1); % +z direction
    % plot the given points
    pastr = "p1 ( " + string(pa(1)) + " , " + string(pa(2)) + " , " + string(pa(3)) + " )";
    text(pa(1)+.5, pa(2)+.5, pa(3)+.5, pastr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(pa(1), pa(2), pa(3), 'ko', 'markersize', 15, 'linewidth', 4);
    pbstr = "p2 ( " + string(pb(1)) + " , " + string(pb(2)) + " , " + string(pb(3)) + " )";
    text(pb(1)+.5, pb(2)+.5, pb(3)+.5, pbstr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(pb(1), pb(2), pb(3), 'ko', 'markersize', 15,'linewidth', 4);
    
    pcstr = "p3 ( " + string(pc(1)) + " , " + string(pc(2)) + " , " + string(pc(3)) + " )";
    text(pc(1)+.5, pc(2)+.5, pc(3)+.5, pcstr,'fontweight','bold','fontsize',14, 'color', 'b');
    plot3(pc(1), pc(2), pc(3), 'ko', 'markersize', 15,'linewidth', 4);
    
    %plot the displacement vector pa to pb
    plot3([pa(1), pb(1)], [pa(2), pb(2)], [pa(3), pb(3)], 'g', 'linewidth', 5);
    plot3([pa(1), pc(1)], [pa(2), pc(2)], [pa(3), pc(3)], 'g:', 'linewidth', 3);
    plot3([close(1,1), pc(1)], [close(1,2), pc(2)], [close(1,3), pc(3)], 'r:', 'linewidth', 2);
    %plot3([close(2,1), pc(1)], [close(2,2), pc(2)], [close(2,3), pc(3)], 'r:', 'linewidth', 2);
end

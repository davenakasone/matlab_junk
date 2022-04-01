%{
    ch14.1 line integration in complex plane
    
    #'a'    ex5     must use integration by path...same start and stop point
    #'b'    ex6     circle centered anywhere
    #'c'    ex7     end points and path matter ...different path = different integral
    #'d'    ex8     ML test    pay attention to the Z and temp Z use
    #1      p1    
    #2      p2  
    #3      p3
    #4      p4
    #5      p5
    #6      p6
    #7      p7
    #8      p8
    #9      p9
    #10     p10
    #11     p11
    #12     p12
    #13     p13   parameterizing a circle from any center
    #14     p14   
    #21     p21
    #22     p22
    #23     p23
    #24     p24
    #25     p25
    #26     p26
    #27     p27
    #28     p28
    #29     p29
    #30     p30


%}
clc;
close all;
clearvars;


                sel = 30;  % CHANGE CHANGE CHANGE
  
                
global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Z; syms Z; Z = X + 1j*Y;% just use it as a hold and imply  when needed
global Zt; syms Zt;            % your temporary Z for integration and differentiation
global Ur; syms Ur; assume(Ur, 'real');   % Ur ( X, Y)  ...real part of Z = Ur(X,Y) + j Vi(X,Y)      
global Vi; syms Vi; assume(Vi, 'real');   % Vi ( X, Y)  ...imag part of Z = Ur(X,Y) + j Vi(X,Y)     
% your function = Ur + j Vi   OR your function = operation(Z)
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve
global pu; syms pu; assume(pu, 'real'); % paramater u for surface trace
global pv; syms pv; assume(pv, 'real'); % paramater v for surface trace

%------------------------------------------------------------------------------------------ #'a'
if sel == 'a'   % one of the most important line integrals
    % int( f(z) = 1/z , 0 2pi)     around unit circle
    syms z;
    path = cos(pt) + 1j * sin(pt);  % parameterize the path around unit circle
        fprintf('path is defined by: %s \n', path);                         % aka exp(jt)
    path_d = diff(path, pt, 1);      % get the derivative of the path
        fprintf('derivative of parameterized path: %s\n', path_d);          % aka  j exp(jt)
    fun = 1/z;
    fun = subs(fun, z, path); % need to substitute path into function to integrate
        fprintf('function in terms of parameter: %s\n', fun);                % aka exp(-jt)
    intg = fun * path_d;  % integrand is parameterized function and derivative of path = 1j
    % because  exp(-jt) { j exp(jt) } = 1j
    % not just integrate
    area = int(intg, pt, 0, 2*sym(pi));  % start on 0, end 2pi...ccw direction
    fprintf('\nthe integral is: %s  \n', area); % j 2pi
    % also note same path, reverse direction is negative:
    area = int(intg, pt, 2*sym(pi), 0);  % start on 0, end 2pi...ccw direction
    fprintf('the integral reversed is: %s  \n', area); % -j 2pi
end


%------------------------------------------------------------------------------------------ #'b'
if sel == 'b'
    z0 = 3 + 1j*3;  % desired center
    rho = abs(z0);  % implied radius
    path = z0 + rho * ( cos(p_t) + 1j*sin(p_t) );  % path of integration  "  z0 + rho exp(jt) "
    path_d = diff(path, p_t, 1); % derivative of path
    
    m = -1; % CHANGE  ...if not -1, area is 0
    syms z;
    fun = (z - z0)^m;
    fun_pt = subs(fun, z, path);  % substitute parameterized path into funciton
    intg = fun_pt * path_d;       % establish integrand
    
    area = int(intg, p_t, 0, 2*sym(pi));
    fprintf('the integral produced: %s\n', area);  %  j 2pi  if m = -1, 0 otherwise  
end


%------------------------------------------------------------------------------------------ #'c'
if sel == 'c'
    syms z;
    fun = real(z);
    
    C = p_t + 1j*2*p_t;  % hypotonuse [0,1]
    C_d = diff(C, p_t, 1); % differentiate this path
    intg = subs(fun, z, C) * C_d;  % obtain the integrand
    area1 = int(intg, p_t, 0, 1);
    fprintf('this line integral is: %s\n', area1);
    
    C1 = p_t;
    C1_d = diff(C1, p_t, 1);
    intgA = subs(fun, z, C1) * C1_d;
    areaA = int(intgA, 0, 1);
    
    C2 = 1 + 1j*p_t;
    C2_d = diff(C2, p_t, 1);
    intgB = subs(fun, z, C2) * C2_d;
    areaB = int(intgB, p_t, 0, 2);
    fprintf(' this other path produces: %s\n', areaA + areaB);  % different becuase of path
end


%------------------------------------------------------------------------------------------ #'d'
if sel == 'd'
    fun = Z^2;
    reslt = fun_rCRtest(fun); % even though it is analytic, still going to find bounds
    integ = int( subs(fun,Z, Zt), Zt);
    integ = subs(integ, Zt, Z);
    ll = subs( integ, Z, 0 );   % can't do subs(integ, [X,Y] , [0,0]
    ul = subs( integ, Z, 1+1j );
    ftc = ul - ll;
    fprintf('\nofficial answer: %s\n', ftc);
    
    pathStart = 0;
    pathStop = 1 + 1j;
    pathL = sqrt(( real(pathStart) - real(pathStop) )^2 + ( imag(pathStart) - imag(pathStop))^2 );
    %fprintf('\nthe path length: %.3f \n', pathL);
    maxM = 2; % abs(fun) can be no greater than 2 on this path
    ML = maxM * pathL; % bound established
    intAbs = abs(ftc);
    fprintf('    the ML bound is %.3f for this path, the abs of integral is: %.3f\n', ML, intAbs);
    % range is good
end


%------------------------------------------------------------------------------------------ #1
if sel == 1
    % z(t) = x(t) + j y(t)  __  t + 1/2 tj  2<t<5
    parmz = [ pt , (1/2)*pt,0];
    rng = [2,5];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(31,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    % z(t) =   t + j(1-t)  0<t<3
    eqn = 3 + 1j + (1-1j)*pt
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,3];
        %fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % z(t) =   t + jt^2  1<t<2
    eqn = pt + 1j*2*pt^2;
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [1,2];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    % z(t) =   t + jt^2  -1<t<1
    eqn = pt + 1j*(1-pt)^2;
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [-1,1];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end

%------------------------------------------------------------------------------------------ #5
if sel == 5
    % z(t) =   3-j+sqrt(10)*exp(-jt)  0<t<2pi
    eqn = 3-1j+sqrt(10)*exp(-1j*pt);
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,2*pi];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % z(t) =   1 + j + exp(-j pi t)  0<t<2
    eqn = 1 + 1j + exp( -pt * 1j * pi);
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,2];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    % z(t) =   2 + 4*exp( pi j t / 2)  0<t<2
    eqn = 2 + 4*exp( pt * 1j * pi / 2);
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,2];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % z(t) =   2 + 4*exp( pi j t / 2)  0<t<pi/2
    eqn = 5*exp( -pt * 1j);
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,pi/2];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end

%------------------------------------------------------------------------------------------ #9
if sel == 9
    % z(t) =   t + j t^3  -2<t<2
    eqn = pt + 1j*pt^3;
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [-2,2];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    % z(t) =   2cos(t) + j sin(t)      0 < t < 2 pi
    eqn = 2*cos(pt) + 1j*sin(pt);
    parmz = [ real(eqn) , imag(eqn), 0];
    rng = [0,2*pi];
        fun_graph_spaceC(parmz, rng);
    pts = zeros(4,1);
    ct = 1;
    for k = rng(1):.1:rng(2)
        temp = subs(parmz(1), pt, k) + 1j*subs(parmz(2), pt, k);
        pts(ct) = temp;
        ct = ct + 1;
    end
        fun_graph2Dcomplex(pts);
        display(pts);
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    % segment (-1,1) --> (1,3)   just use eqn of line y = mx + b
    ptA = [-1,1];
    ptB = [1,3];
    slope = ( ptB(2) - ptA(2) ) / ( ptB(1) - ptA(1) );
    intp = ptA(2) - slope*ptA(1);
    syms x;
    syms y;
    y = slope*x + intp;
    fprintf(' use  y =  %s\n', y);
    fprintf(' or x = t, y = j(t + 2),    -1<t<1 ...natural\n');
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    % (0,0) --> (2,1)  along axes....2 seperate ones
    fprintf('first use  x = t , y = 0j   0<t<2\n');
    fprintf('then use x = 2, y = t   0<t<1\n');
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    % upper half abs( z - 2 + j) = 2     (4,-1) to (0, -1)
    % ...  (x-2)^2 + (y+1)^2 = 4    x-2 = 2 cos(t)    y+1 = 2 sin(t)
    fprintf('path is 2 - j + 2*exp(jt)\n');
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    fprintf( 'unit circle cw is   cos(t) + j sin(t)     2pi --> 0\n');
    fprintf( '    or use   cos(t) - j  sin(t)  0,2pi\n');
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    %  x^2 -4y^2 = 4   branch through (2,0)
    %  x^2 / 4  - y^2 = 1   --> x^2 / 4 - 1 = y^2    , y = sqrt(  x^2 \4  - 1)
    % could go natural, or x = 2cost, y = sint
    % know how to reference all the parameters for 
end
    

%------------------------------------------------------------------------------------------ #21
if sel == 21
    % integrating real(z) (1+j) to (3+3j)
    fun = real(Z);
        %fun_rCRtest(fun); % FAILED...can't use ftc method
    path = pt + 1j*pt; % on the y = x line  [ 1, 3 ]
    path_d = diff(path, pt, 1);  % 1 + j
    fun = real(path); %  f(z(t) = pt   after sub
    intg = fun * path_d; %  pt + j pt
    area = int(intg, pt, 1, 3);
    fprintf(' integral is:  %s\n', area);  %  4 + 4i
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    % integrating real(z) (1+j) to (3+3j)   on parabola y = 1 + 1/2(x-1)^2   go natural
    fun = real(Z); % already know it is not analytic
    path = pt + 1j*( 1 + .5*(pt-1)^2);  
    path_d = diff(path, pt, 1); %   1 + j(t-1)
    fun = real(path); % sub  f(z(t)) = pt
    intg = fun * path_d; % = t + j (t^2 -t)
    area = int(intg, pt, 1, 3);
    fprintf('integral is :  %s \n', area);   % 4 + 14i/3 
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    % exp(z)  pij to 2pij
    fun = exp(Z);
        fun_rCRtest(fun);  % it is analytic...just use ftc
     intg = exp(pt);
     area = exp(2*1j*pi) - exp(1j*pi)
     %fprintf('\nintegral is: %d \n', area);
     
     % parameterize just to check...it is a mess
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    % cos(2z) on semi circle abs(z) = pi  ... x>0    -pij to pij
    fun = cos(2*Z);
    %fun_rCRtest(fun);  % it's analytic, use ftc
    fun = cos(2*Zt);
    antiD = int(fun, Zt)
    ll = subs(antiD, Zt, -1j*pi);
    ul = subs(antiD, Zt, 1j*pi);
    area = ul - ll; %     aka  " j sinh(2pi) "
end


%------------------------------------------------------------------------------------------ #25
if sel == 25
    % z exp(z^2)    1 on axis to j    (1,0) --> (0,j)
    fun = Z * exp(Z^2);
        fun_rCRtest(fun); % it is analytic, but an abortion
        %fun_graphSurCon(fun, 5); %  it is an abortion
    
    antiD = (1/2)*exp(pt^2);
    ll = subs(antiD, pt, 1);
    ul = subs(antiD, pt, 1*j);
    area = simplify(ul - ll)
        
    path = pt + 1j*(1-pt); % path is 1<t<0
    path_d = diff(path, pt, 1);  % 1 - j
    fun = path * exp(path^2); % sub
    intg = fun * path_d;
    area = simplify(int(intg, pt, 1, 0))     %  - exp(-1)/2 + exp(1)/2   "-sinh(1)"
    double(simplify(area + sinh(1)))
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    fun = Z + Z^(-1);   % 1/z should give it away it is not analytic in unit circle
    fun_rCRtest(fun); % not analytic
    %fun_graphSurCon(fun, 5);  % that is a flithy looking function
    
    path = cos(pt) + 1j*sin(pt); % ccw, 0,2pi  ...exp(jt)
    path_d = diff(path, pt, 1);  %  -sin(t) + j cos(t)  --> j exp(jt)
    subF = path + 1/path;  %   exp(jt) + exp(-jt) ...  2cosh(jt)
    intg = subF * path_d;  %  aka    j exp(2jt) + j
    area = int(intg, pt, 0, 2*sym(pi))   %  pi*2i
end


%------------------------------------------------------------------------------------------ #27
if sel == 27
    % sec(z)^2 is not analytic in some places, but is in this domain
    fun = (sec(Z))^2;
    fun_rCRtest(fun);
    antiD = tan(pt);
    ll = subs(antiD, pt, sym(pi)/4);
    ul = subs(antiD, pt, 1j*sym(pi)/4);
    area1 = simplify(ul-ll)  % tanh(pi/4)*1i - 1
    
    path = pt*sym(pi)/4 + 1j*(1-pt)*sym(pi)/4; % [ 1 , 0]
    path_d = diff(path, pt, 1);
    subF = (sec(path))^2;
    intg = subF * path_d;
    area2 = simplify(int(intg, pt, 1, 0))
    double(area1-area2)
end


%------------------------------------------------------------------------------------------ #28
if sel == 28
    % you can almost tell this one won't be analytic
    fun = (5/(Z-2j))-(6/(Z-2j)^2);
    fun_rCRtest(fun); % it actually passed... but it is on a loop, so no ftc  maybe split?
    %fun_graphSurCon(fun, 5);
    % had that general root property also, to cancel second integral
    
    path = 4*cos(pt) + 1j*( 2 + 4*sin(pt));
    path_d = diff(path, pt, 1);
    subF = (5/(path-2j))-(6/(path-2j)^2);
    intg = subF * path_d;
    area = int(intg, pt, 2*sym(pi),0) 
end


%------------------------------------------------------------------------------------------ #29
if sel == 29
    % f(z) = imag(z^2)   (0,0) -->  (1,0) --> (0,j) -->
    fun = imag(Z^2);
    %fun_rCRtest(fun); % not analytic...going to have to make a line integral for triangle
    
    % (0,0) --> (1,0)
    path = pt + 1j*0;
    path_d = diff(path, pt, 1); %  = 1
    subF = imag(path^2); % should be 0
    intg = subF * path_d;
    seg1 = int(intg, pt, 0, 1);
    
    % (1,0) --> (0,j)
    path = pt + 1j*(1-pt);
    path_d = diff(path, pt, 1); %  = 1 - j
    subF = imag(path^2); 
    intg = subF * path_d;
    seg2 = int(intg, pt, 0, 1);
    
    % (0,j) --> (0,0)
    path = 0 + 1j*(pt);
    path_d = diff(path, pt, 1); %  =  j
    subF = imag(path^2); 
    intg = subF * path_d;
    seg3 = int(intg, pt, 1, 0);
    (seg1+seg2+seg3)*-1  % its cw not ccw  
end


%------------------------------------------------------------------------------------------ #30
if sel == 30
    % f(z) = real(z^2)   (0,0) -->  (0,j) --> (1,j) --> (1,0) --> (0,0)
    fun = real(Z^2);
    %fun_rCRtest(fun); % failed
    
    %(0,0) --> (0,j)
    path = 0 + 1j*pt;
    path_d = diff(path, pt, 1);
    subF = real(path^2);
    intg = subF * path_d;
    seg1 = int(intg, pt, 0, 1);
    
    %(0,j) --> (1,j)
    path = pt + 1j*1;
    path_d = diff(path, pt, 1);
    subF = real(path^2);
    intg = subF * path_d;
    seg2 = int(intg, pt, 0, 1);
    
    %(1,j) --> (1,0)
    path = 1 + 1j*pt;
    path_d = diff(path, pt, 1);
    subF = real(path^2);
    intg = subF * path_d;
    seg3 = int(intg, pt, 1, 0);
    
    %(1,0) --> (0,0)
    path = pt + 1j*0;
    path_d = diff(path, pt, 1);
    subF = real(path^2);
    intg = subF * path_d;
    seg4 = int(intg, pt, 1, 0);
    
    seg1+seg2+seg3+seg4   % - 1 - j
end
    
    
    
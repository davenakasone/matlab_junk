%{
    ch 14.3 Cauchy's Integral Formula
    
    #'a'    ex1
    #'b'    ex2
    #'c'    ex3    must include the point for line integral to be valid
    #1      p1 to p4  check these
    #5      p5
    #6      p6
    #7      p7
    #8      p8
    #10     p10
    #11     p11
    #12     p12
    #13     p13
    #14     p14

%}
clc;
close all;
clearvars;


                sel = 14;  % CHANGE CHANGE CHANGE



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
if sel == 'a'
    fprintf(' should just use formula to get    2*pi*j*e^2  , j %.3f\n',...
        exp(2)*pi*2);
    funZt = exp(Zt) / ( Zt - 2); % see it and immidiatley apply formula
    
    path = cos(pt) + 1j*sin(pt); % radius = 1
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line1 = int(intg, pt, 0, 2*sym(pi)) % not valid , z = 2 not enclosed
    
    path = 3*cos(pt) + 3j*sin(pt); % radius = 3
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line2 = int(intg, pt, 0, 2*sym(pi)) % valid, z = 2 is enclosed and matches expectation
end


%------------------------------------------------------------------------------------------ #'b'
if sel == 'b'
    funZ = (Z^3 -6) / (2*Z -1j);
    %fun_graphSurCon(funZ, 5);
    funZt = (Zt^3 -6) / (2*Zt -1j); % hole at z =  1/2 j
    
    path = cos(pt) + 1j*sin(pt); % radius = 1
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, pt, 0, 2*sym(pi)) % valid, z = 1/2 j is enclosed and matches expectation
end


%------------------------------------------------------------------------------------------ #'c'
if sel == 'c'
    funZ = (Z^2 + 1) / (Z^2 -1);
    %fun_graphSurCon(funZ, 5);
    funZt = (Zt^2 + 1) / (Zt^2 -1); % holes at z =  +/- 1
    
    path = 1+cos(pt) + 1j*sin(pt); % radius = 1 , center (1,0)  hole included
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circA = int(intg, pt, 0, 2*sym(pi))   % pi*2i     m=-1
    
    path = .25+cos(pt) + 1j*sin(pt); % radius = 1 , center (.25,0)  hole included
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circB = int(intg, pt, 0, 2*sym(pi))   % pi*2i     m=-1
    
    path = cos(pt)-1 + 1j*sin(pt); % radius = 1 , center (-1,0)  hole included
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circC = int(intg, pt, 0, 2*sym(pi))   % -pi*2i     m=-1
    
    path = cos(pt) + 1j*sin(pt)+1; % radius = 1 , center (0,1)  NO HOLES INCLUDED
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circD = int(intg, pt, 0, 2*sym(pi))   % 0     m=-1
end

%------------------------------------------------------------------------------------------ #1
if sel == 1
    funZ = Z^2 / (Z^2 -1);           %  holes at z = +/- 1
        %fun_graphSurCon(funZ, 5);
        %fun_rCRtest(funZ);          % it is analytic
    funZt = (Zt^2) / (Zt^2 -1); 
    
    path = -1+cos(pt) + 1j*sin(pt); % radius = 1 , center (-1,0)  hole included
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, 0, 2*sym(pi))   % - pi j  
    
    % radius = pi/2 , center (1,j)  hole included
    r = sym(pi)/2;
    path = (1 + cos(pt))*r + 1j*(1 + sin(pt))*r; 
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ2 = int(intg, pt, 0, 2*sym(pi))   %  pi j   use formula   this isn't working check
    
    % radius 1.4, center (0, -1) has hole
    r = 1.4;
    path = (cos(pt))*r + 1j*(-1 + sin(pt))*r;  % one point inside
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ3 = int(intg, pt, 0, 2*sym(pi))   %  pi j
    
    % radius 1.4, center (0, -1) has hole
    r = 1.4;
    path = (cos(pt))*r + 1j*(-1 + sin(pt))*r;  % one point inside
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ3 = int(intg, pt, 0, 2*sym(pi))   %  pi j 
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    funZt = (cos(3*Zt))/(6*Zt);
    expand(funZt)
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, 0, 2*sym(pi))    % take that 1/6 out, pole at z=0  >> (pi*j)/3
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    check = double(-2*sin(2/sym(pi)) + 1j*2*cos(2/sym(pi))); % good, 2j exp(2j/pi)
    funZt = exp(2*Zt)/(Zt*sym(pi)-1j);
    expand(funZt);
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, 0, 2*sym(pi));   
    double(circ)
    double(check)
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
     % pi/8  is good
    funZt = Zt^3 / (2*Zt - 1j);
    expand(funZt);
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, 0, 2*sym(pi))  
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    funZt = ( (Zt^2) * sin(Zt) ) / (4*Zt - 1j);
    expand(funZt);
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, 0, 2*sym(pi)); % wrong
    double(circ)
    double(1j*sym(pi)*sin(1/4)/32) % correct
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    funZ = (Z^3 - 6) / (Z - (1/2)*1j);
    %fun_graphSurCon(funZ,2);
    
    funZ = (sin(Z)) / (Z - (1/2)*pi);
    fun_graphSurCon(funZ,1);
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    funZ = 1 / ( Z^2 + 4 );
    %fun_graphSurCon(funZ,5);
    funZt = 1 / ( Zt^2 + 4 );
    
    rotz = fun_demovFun();
    %rotz    %the 2 holes
    
    % circle, ccw, (0, 2) , r = 2
    centR = 0;
    centI = 2;
    r = 2;
    start = 0;
    stop = 2*sym(pi);
    path = r*(centR + cos(pt)) + 1j*r*(centI + sin(pt));
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ = int(intg, pt, start, stop); % failed
    
    f = 1 / ( Z + 2j);
    z0 = 2j;
    fSub = subs(f, Z, z0);
    circ = fSub * 2j*sym(pi); % correct: pi/2
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    funZ = Z / ( Z^2 + 4*Z +3 );
    %fun_graphSurCon(funZ,5);
    funZt = Zt / ( Zt^2 + 4*Zt +3 );
    
    rotz = fun_demovFun();
    %rotz    %the 2 holes
    
    % circle, ccw, (-1, 0) , r = 2
    centR = -1;
    centI = 0;
    r = 2;
    start = 0;
    stop = 2*sym(pi);
    path = r*(centR + cos(pt)) + 1j*r*(centI + sin(pt));
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop) % worked   pi*2j  correct
    
    f = Z / ( Z + 3);
    z0 = -1;
    fSub = subs(f, Z, z0);
    circ2 = fSub * 2j*sym(pi) % failed: -j pi
    
    % pf--> m = -1 works    can't use formula because there are 2 holes..not 1
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    funZ = ( Z + 2 ) / ( Z - 2 );
    %fun_graphSurCon(funZ,5);
    funZt =  (Zt + 2 ) / ( Zt - 2 );
    
    rotz = fun_demovFun();
    %rotz    %the 2 holes
    
    % circle, ccw, (1, 0) , r = 2
    centR = 1;
    centI = 0;
    r = 2;
    start = 0;
    stop = 2*sym(pi);
    path = r*(centR + cos(pt)) + 1j*r*(centI + sin(pt));
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop) % worked   pi*8i  correct
    
    f = Z + 2;
    z0 = 2;
    fSub = subs(f, Z, z0);
    circ2 = fSub * 2j*sym(pi) % worked  pi*8i
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    funZ = exp(Z) / ( Z*exp(Z) - 2j*Z );
    %fun_graphSurCon(funZ,5);
    funZt =  exp(Zt) / ( Zt*exp(Zt) - 2j*Zt );
    
    %rotz = fun_demovFun();
    %rotz    %the 2 holes         only 0 is inside, ln(2) + jpi/2 are outside
    
    % circle, ccw, (0, 0) , r = .6
    centR = 0;
    centI = 0;
    r = .6;
    start = 0;
    stop = 2*sym(pi);
    path = r*(centR + cos(pt)) + 1j*r*(centI + sin(pt));
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop); % ?
    
    
    f = exp(Z) / (exp(Z) - 2j) ;
    z0 = 0;
    fSub = subs(f, Z, z0);
    circ2 = fSub * 2j*sym(pi); % pi*(- 4/5 + 2i/5)
    simplify(circ2);
end
    



    
%{
    ch 14.2 Cauchy's Integral Theorem
    
    #1  p1    theorem holds for any closed path
    #2  p2    connectedness
    #6  p6    path independence
    #7  p7    deformations, still integrable if pole outside 
    #8  p8    partial fractions are allowed
    #9  p9    to use or not to use Cauchy's integral theorem
    #10 p10     "  "
    #11 p11     "  "
    #12 p12     "  "
    #13 p13     "  "
    #14 p14     "  "
    #15 p15     "  "
    #16 p16     "  "
    #17 p17     "  "
    #18 p18     "  "
    #19                 taylor series
    #20
    #21
    #22
    #23

%}
clc;
close all;
clearvars;


                sel = 23;  % CHANGE CHANGE CHANGE



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


%------------------------------------------------------------------------------------------ #1
if sel == 1
    % verify  f(z) = z^2 over square  +/- 1 , +/- j vertecies is 0
    fun = Z^2;
    %fun_rCRtest(fun);  % it is analytic and there is nothing to suggest it is undefined
    fun = Zt^2;
    
    seg1a = int(fun, Zt, 1 - 1j, 1 + 1j);
    seg2a = int(fun, Zt, 1 + 1j, -1 + 1j);
    seg3a = int(fun, Zt, -1 + 1j, -1 - 1j);
    seg4a = int(fun, Zt, -1 - 1j, 1 - 1j);
    circA = seg1a+seg2a+seg3a+seg4a
    
    seg1b = int(fun, Zt, 1 - 1j, -1 - 1j);
    seg2b = int(fun, Zt, -1 - 1j, -1 + 1j);
    seg3b = int(fun, Zt, -1 + 1j, 1 + 1j);
    seg4b = int(fun, Zt, 1 + 1j, 1 - 1j);
    circB = seg1b+seg2b+seg3b+seg4b
    
    % now try a circle              any path will be 0
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = path^2;
    intg = subF * path_d;
    circC = int(intg, pt, 0, 2*sym(pi))
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    % f(z) = 1/z is going to be a problem if z = 0, aslo, antiD = Ln(z) so any neg real
    % try (1,j) --> (3, j) --> (2, 2j) --> (1,j)  or anything without z=0
    funZ = 1/Zt;
    seg1 = int(funZ, Zt, 1 + 1j, 3 + 1j);
    seg2 = int(funZ, Zt, 3 + 1j, 2 + 2j);
    seg3 = int(funZ, Zt, 2 + 2j, 1 + 1j);
    circ = double(simplify(seg1+seg2+seg3))  % almost 0
    
    % f(z) = exp(1/z^2) / ( z^2 + 16 )   good for all except z = 0, 4j, -4j
    
    % can always split into segments to get around the poles
    % ie  f(z) = 1/z^4 + 1   has 4 poles, so it is 5-connected
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % f(z) = exp(z)    (0,0) --> (1,j)     vs    (0,0) --> (1,0) --> (1,j)
    funZ = exp(Zt);
    
    test1 = int(funZ, Zt, 0, 1 + 1j)           % same start/stop = same integral, any path
    
    seg1 = int(funZ, Zt, 0, 1);
    seg2 = int(funZ, Zt, 1, 1 + 1j);
    test2 = seg1+seg2
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    %  f(z) = 1 / ( z^2 + 4)  on abs(z-2) = 2
    % poles are z = 2j , -2j      abs( 2j -2) and abs(-2j-2) = 2 sqrt(2)  ...outside
    % since poles outside domain of interest, it will work:
    funZ = 1 / (Zt^2 + 4 ); 
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = 1 / (path^2 + 4 );
    intg = subF * path_d;
    circ = int(intg, pt, 0, 2*sym(pi))
    % not happening if domain gets too big and contains poles
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % f(z) = ( 2z + j3) / (z^2 + 1/4) --> poles at .5
    funZ = ( 2*Zt + 3j) / ( Zt^2 + (1/4));
    rad = .2; % CHANGE ...keep under .5
    
    path = rad*cos(pt) + rad*1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Zt, path);
    intg = subF * path_d;
    circ = int(intg, 0, 2*sym(pi))
    
    % f(z) = (z+1)/(z^2 + 2z)   poles at -2, 0
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    funZ = exp(-Z^2);
    %fun_rCRtest(funZ);  % it is analytic, should be 0 around unit circle
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi)) % it is 0
    % by ftc not possible, same start and stop
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    funZ = tan( (1/4) * Z);
    %fun_rCRtest(funZ); % it is analytic, but it is not defined at some points outside
    % should be 0
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi))  % it is 0, as expected
end

%------------------------------------------------------------------------------------------ #11
if sel == 11
    funZ = 1 / ( 2*Z - 1);
    %fun_rCRtest(funZ); % not analytic, and has a pole at Z = 1/2
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi))  % it is 0 but this is wrong
    % have to use (z-z0)^m   = 2pi j if m = -1   
    % should not have even tried integrating...use formula
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    funZ = conj(Z)^3;
    %fun_rCRtest(funZ); % not analytic
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi))  % says 0
    expand(conj(X+1j*Y)^3); % no holes  X^3 - X^2*Y*3j - 3*X*Y^2 + Y^3*1j
    % use conj(z) = exp(-jt)  ....conj(z)^3  = exp(-3jt)
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    funZ = 1 / ( Z^4 - 1.1);
    %fun_rCRtest(funZ); % analytic, but will have some holes around z^4 = 1.1  4 roots
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi));  % says 0 , good since roots outside
    x=fun_demovFun();  % all more than 1
    x(1) % 1.0241
    x(2) % 0.0000 + 1.0241i
    x(3) % -1.0241
    x(4) % 0.0000 - 1.0241i
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    funZ = 1 / conj(Z);
    %fun_rCRtest(funZ); % not analytic
    % but it  hole at Z = 0 and it does have a -1 power...so it should be 2pi j
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZ, Z, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi))  % 0, by multiply connected domains
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    funZ = imag(Z);
    %fun_rCRtest(funZ); % not analytic
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    fun = imag(Zt);
    subF = subs(fun, Zt, path);
    intg = subF * path_d;
    lineInt = int(intg, pt, 0, 2*sym(pi))  % -pi
end

%------------------------------------------------------------------------------------------ #16
if sel == 16
    funZ = 1 / (sym(pi)*Z - 1);
    %fun_rCRtest(funZ); %  there is a clear hole at Z = 1/pi
    fun = 1 / (sym(pi)*Zt - 1);
    % you have interation across regions or m = -1
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(fun, Zt, path);
    intg = subF * path_d;
    lineINt = int(intg, pt, 0, 2*sym(pi)) %  2j
    % or since m = -1 and area is 2 pi j   with 1/pi in front, 2j confirmed
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    funZ = 1 / abs(Z)^2; 
    %fun_rCRtest(funZ); % not analytic
    funZt = 1 / abs(Zt)^2;   % clear hole at Z = 0
    
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line = int(intg, pt, 0, 2*sym(pi))  % 0, by connected domain
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    funZ = 1 / (4*Z - 3);   % hole at Z = 3/4...1/4 outside, should be 2pi j / 4
    %fun_rCRtest(funZ); % not analytic anyway
    funZt = 1 / (4*Zt - 3);
    
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line = int(intg, pt, 0, 2*sym(pi))  % (pi*1i)/2      works
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    funZ = cot(Z)*Z^3;     % holes where sin(Z) = 0, but outside, except Z = 0
    %fun_rCRtest(funZ); % fail
    funZt = cot(Zt)*Zt^3;
    
    path = cos(pt) + 1j*sin(pt);   % had to use taylor series to see if analytic  wtf
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line = double(simplify(int(intg, pt, 0, 2*sym(pi))))  %almost 0
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    funZ = log(1-Z);     % holes if Z = 1, but not in boundary   or any 1-Z < 0 real
    %fun_rCRtest(funZ); % fail
    funZt = log(1-Zt);
    
    path = pt + 1j*(2*pt-1);   % (0,-j) to (1,j)
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line1 = int(intg, pt, 0, 1);  % [0,1]
    
    path = pt + 1j;   % (1,j) to (0,j)
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line2 = int(intg, pt, 1, 0); % [1,0]
    
    path = pt + 1j*(2*pt+1);   % (0,j) to (-1,-j)
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line3 = int(intg, pt, 0, -1); % [0,-1]
    
    path = pt + -1j;   % (-1,-j) to (0,-j)
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line4 = int(intg, pt, -1, 0); % [-1,0]
    
    double(simplify(line1+line2+line3+line4)) % call it 0
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    funZ = 1 / (Z-3j);
    %fun_rCRtest(funZ); % fails...but m = -1     going to be 2pi j
    funZt = 1 / (Zt-3j);
    
    path = sym(pi)*cos(pt) + sym(pi)*1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line1 = int(intg, pt, 0, 2*sym(pi)) % it is  
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    funZ = real(Z);
    funZt = real(Zt);
    %fun_rCRtest(funZ); % failed
    
    path = cos(pt) + 1j*sin(pt);
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line1 = int(intg, pt, 0, sym(pi));
    
    path = pt;
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    line2 = int(intg, pt, -1, 1);
    
    line1 + line2
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    funZ = (2*Z - 1) / (Z^2 - Z);
    funZt = (2*Zt - 1) / (Zt^2 - Zt);
    fun_rCRtest(funZ); % technicall passes, but has holes at Z = 1, 0
    % pf -->  1/z-1 + 1/z   m = -1    4pi j   ( 2 * 2pi j)    #24 is the same...see m=-1
end
    
    
    
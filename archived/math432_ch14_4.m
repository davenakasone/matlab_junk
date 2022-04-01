%{
    ch 14.4 Derivatives of Analytic Functions
    
    #1  
    #2
    #4
    #5
    

%}
clc;
close all;
clearvars;


                sel = 5;  % CHANGE CHANGE CHANGE



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
    funZt = sin(Zt) / Zt^4;  % notice pole at z = 0
    
    % taking it ccw on r = 1, center (0,0) 'unit circle'
    centR = 0;
    centI = 0;
    r = 1;
    start = 0;
    stop = 2*sym(pi);
    path = r * ( centR + cos(pt) ) + 1j * r * ( centI + sin(pt) );
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop)  % -(pi*1i)/3   as expected , hole contained
    
    order = 3; % 3rd derivative will solve this
    z0 = 0;    % hole at z = 0 implied
    mult = sym(pi)*j*2 / factorial(order);
    f = sin(Zt); % implied f(z)
    dfn = diff(f, Zt, order);
    fz0 = subs(dfn, Zt, z0);
    circ2 = mult * fz0               %  -(pi*1i)/3   as expected
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    funZt = Zt^6 / (2*Z-1)^6;  % notice pole at z = 1/2  
    
    % taking it ccw on r = 1, center (0,0) 'unit circle'
    centR = 0;
    centI = 0;
    r = 1;
    start = 0;
    stop = 2*sym(pi);
    path = r * ( centR + cos(pt) ) + 1j * r * ( centI + sin(pt) );
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop)  % fails, holes?
    
    order = 5; % 5th derivative will solve this
    z0 = 1/2;    % hole at z = 1/2 implied
    mult = sym(pi)*j*2 / factorial(order);
    f = (1/2^6)*(Zt^6); % implied f(z)
    dfn = diff(f, Zt, order);
    fz0 = subs(dfn, Zt, z0);
    circ2 = mult * fz0               %  (pi*3i)/32   as expected   PF a lot better ....m=-1
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    funZt = exp(Zt)*cos(Zt)/(Zt-sym(pi)/4)^3;  % notice pole at z = pi/4  
    
    % taking it ccw on r = 1, center (0,0) 'unit circle'
    centR = 0;
    centI = 0;
    r = 1;
    start = 0;
    stop = 2*sym(pi);
    path = r * ( centR + cos(pt) ) + 1j * r * ( centI + sin(pt) );
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop)  % probably good
    
    order = 2; % 2nd derivative will solve this
    z0 = sym(pi)/4;    % hole at z = pi/4 implied
    mult = sym(pi)*j*2 / factorial(order);
    f = exp(Zt)*cos(Zt); % implied f(z)
    dfn = diff(f, Zt, order);
    fz0 = subs(dfn, Zt, z0);
    circ2 = mult * fz0               %  -2^(1/2)*pi*exp(pi/4)*1i   as expected   PF a lot better ....m=-1
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    funZt = cosh(2*Zt) /(Zt-(1/2))^4;  % notice pole at z = 1/2  
    
    % taking it ccw on r = 1, center (0,0) 'unit circle'
    centR = 0;
    centI = 0;
    r = 1;
    start = 0;
    stop = 2*sym(pi);
    path = r * ( centR + cos(pt) ) + 1j * r * ( centI + sin(pt) );
    path_d = diff(path, pt, 1);
    subF = subs(funZt, Zt, path);
    intg = subF * path_d;
    circ1 = int(intg, pt, start, stop)  % probably good
    
    order = 3; % 3rd derivative will solve this
    z0 = 1/2;    % hole at z = 1/2 implied
    mult = sym(pi)*j*2 / factorial(order);
    f = cosh(2*Zt); % implied f(z)
    dfn = diff(f, Zt, order);
    fz0 = subs(dfn, Zt, z0);
    circ2 = mult * fz0               %  (pi*sinh(1)*8i)/3   as expected  
end
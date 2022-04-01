%{

    Chapter 15.4 Taylor and Maclaurin Series     ...they are power series

    #'a'    basic demo...use inf expasion point, Laurent series comes out  even multi var possible
    #3      taylor  sin(2z^2)
    #4      taylor z+2 / 1-z^2
    #5      taylor 1 / 2+z^4       equate to known 1/1-z
    #6      1/1+3jz
    #7      (cos(1/2)z)^2   reduce power
    #8      (sinZ)^2        ""
    #9      starts with integral...go term by term
%}
clc;
close all;
clearvars;
sympref('MatrixWithSquareBrackets', 1);
sympref('PolynomialDisplayStyle', 'ascend');
    old_val = sympref('HeavisideAtOrigin', 1);


                    sel = 9;  % CHANGE CHANGE CHANGE

global Z; syms Z            % your temporary Z for integration and differentiation
global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Zxy; syms Zxy; Zxy = X + 1j*Y;% just use it as a hold and imply  when needed
global in; syms in; assume(in, {'real', 'integer'}); % index n
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
%------------------------------------------------------------------------------------------ #'a'
if sel == 'a'
    ctr = 0;     % it is MacClurin
    fun = sin(Z);
    terms = 7;
    
    % savage way 1
    syms sers;
    syms coef;
    sers = subs(fun, Z, ctr);    % take f(z0)  
    for n = 1:terms
        coef = ( (Z-ctr)^n / factorial(n) ) * subs(diff(fun, Z, n), Z, ctr);
        sers = sers + coef;
    end
        %pretty(sers);
    sers = subs(sers, Z, Zxy);
    fun = subs(fun, Z, Zxy);
        %fun_graphSurConX([fun,sers], 1);
    delete(sers);
    
    %savage way 2
    sers = symsum(     ( ((-1)^in) * ( (Z^(2*in+1)) / factorial(2*in+1)) )  , in, 0, terms);
        %pretty(sers);
    sers = subs(sers, Z, Zxy);
        %fun_graphSurConX([fun,sers], 1);
    delete(sers);
    
    % best
    fun = sin(Z);
        %pretty(taylor(fun));
        %pretty(taylor(fun, 'ExpansionPoint', 2)); % center  /  z0
        %pretty(taylor(fun, Z, 'ExpansionPoint', 2));  % give it a variable explicitly
        %t10 = taylor(fun, Z, 'Order', 10, 'ExpansionPoint', 0)
        %t = taylor(fun, Z, 'Order', 90, 'ExpansionPoint', 0, 'OrderMode', 'relative') % for accuracy
    sers = taylor(fun, Z, 'order', 13, 'ExpansionPoint', 0);
    sers = subs(sers, Z, Zxy);
    fun = subs(fun, Z, Zxy);
    fun_graphSurConX([fun, sers], 1);   
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    fun = sin(2*Z^2);
    sers = taylor(fun, Z, 'Order', 15, 'ExpansionPoint', 0);
    
    seq = ((-1)^in)*(2^(2*in+1))/factorial(2*in+1);
    an = seq;
    an1 = subs(seq, in, in+1);
    temp = limit( abs(an/an1), in, inf);
    fprintf('looks like it converges on any radius\n');
    actual = symsum(seq*Z^(4*in+2), in, 0, inf);
    
    sers = subs(sers, Z, Zxy);
    fun = subs(fun, Z, Zxy);
        %fun_graphSurConX([fun, sers], 1);
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    fun = (Z+2)/(1-Z^2);
    funZ = subs(fun, Z, Zxy);
        %fun_graphSurConX(funZ, 3);
    t = taylor(fun, Z, 'Order', 20, 'ExpansionPoint', 0);
    t = subs(t, Z, Zxy);
        %fun_graphSurConX([funZ,t], 1);
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    fun = 1/(2+Z^4);
    funZ = subs(fun, Z, Zxy);
        %fun_graphSurConX(funZ,5);
    t = taylor(fun, Z, 'Order', 7, 'ExpansionPoint', 0);
    % set the form , factor out 1/2 and use f(z) = 1/1-z as known series
    seq = (-1)^in * (1/2)^(in+1);
    sers = seq * Z^(4*in);
    an = seq;
    an1 = subs(seq, in, in+1);
    R = limit( abs(an/an1), in, inf);
    R = R^(1/4); % take the 1/4 root because Z^4n
    sum = symsum(sers, in, 0, inf);
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    fun = 1/(1+(3j*Z));
    % no need to factor, just equate to known 1/1-z  ...z = -3jZ here
    an = (-3j)^in;
    an1 = subs(an, in, in+1);
    R = limit( abs(an/an1), in, inf)
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    fun = (cos((1/2)*Z))^2;
    fun = (1/2) + (1/2)*cos(Z); % helps to reduce
    an = (-1)^in * (1/factorial(2*in));
    an1 = subs(an, in, in+1);
    R = limit ( abs(an/an1), in, inf) % oo as expected....cos good for any
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    fun = (sin(Z))^2;
    fun = (1/2)-(1/2)*cos(2*Z); % reduced
    an = (1/2)*(-1)^in * (1/(factorial(2*in+1)));
    an1 = subs(an, in, in+1);
    R = limit( abs(an/an1), in, inf); % oo as expected
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
end
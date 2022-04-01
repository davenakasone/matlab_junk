%{
    chapter 15.3  Functions from power series

    #'a' ex1
    #3   see the trick  ln(f(z))
    #5   find R
    #6   "        "   need some algebra to set coeff z^n       and z^2n
    #7   sneaky z^2
    #8   " "
    #9   " "

%}
clc;
close all;
clearvars;
        
                      sel = 9;  % CHANGE CHANGE CHANGE



global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Z; syms Z; %Z = X + 1j*Y;% just use it as a hold and imply  when needed
global Zt; syms Zt;            % your temporary Z for integration and differentiation
global Ur; syms Ur; assume(Ur, 'real');   % Ur ( X, Y)  ...real part of Z = Ur(X,Y) + j Vi(X,Y)      
global Vi; syms Vi; assume(Vi, 'real');   % Vi ( X, Y)  ...imag part of Z = Ur(X,Y) + j Vi(X,Y)
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
global T; syms T; assume(T, 'real');  % period bound, continous case
global in; syms in; assume(in, {'real', 'integer'}); % index n
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect

%------------------------------------------------------------------------------------------ #'a'
if sel == 'a'
    syms a;
    %a = (1,N); % coeff vector
    seq = Z^in / factorial(in);
    sers = symsum(seq, in, 0, N);
    temp = subs(sers, [Z,N], [3,9])-exp(3);
    double(temp)
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    f = in^(1/in);
    ln_f = log(in)/in;    %   now it is  ln(f(n)) = ln(n)/n  
    l = limit(ln_f, in, inf);  % it is 0, but that means if ln(f(n)) = 0    lim = exp0 = 1
    lm = exp(l)
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    coeff = in*(in-1) / (2^in);
    var = (Zt-2j)^in;
    seq = coeff * var;
    ser = symsum(seq, in, 2, N);
    
    an = coeff;
    an1 = subs(an, in, in+1);
    R = limit( an/an1, in, inf)
    
    % any derivative will work also and be R = 2 ...just bump the index
    d_an = in*coeff;
    d_an1 = subs( d_an, in, in+1);
    Rd = limit( d_an/d_an1, in, inf)
    
    
    %{
    don't chop up the coeff and z^n term... must isolate to (z - z0)^n form
    seq = ( in*(in-1)/(2^in) ) * (Zt - 2j)^in;
    ser = symsum(seq, in, 2, N);
        %lim = simplify(limit(ser, N, inf))   % not too helpful
    %  Cauchy  Hadamard:  R = 1 / L* = lim n->oo abs( an / an+1)
    an = in*(in-1);  % leave z^n  =  ( z/2  - j)^n
    an1 = subs(an, in, in+1);
    R = limit( an/an1, in, inf)  % R = 1
    %} 
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    temp1 = ((-1)^(in));
    temp2 = (1/(2*in+1));
    temp3 = ((1/(2*sym(pi)))^(2*in+1));
    an = temp1*temp2*temp3;
    an1 = subs(an, in, in+1);
    R_lim = limit( abs(an/an1), in, inf);  %  4*pi^2
    R = sqrt(R_lim); %  2*pi
    
    d_an = ((-1)^in) / ( 2*sym(pi))^(2*in+1);
    d_an1 = subs(d_an, in, in+1);
    Rd_lim = limit( abs(d_an/d_an1), in, inf); % 4*pi^2
    Rd = sqrt(Rd_lim);  %  2*pi
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    an = in/(3^in);
    an1 = subs(an, in, in+1);
    R = limit( abs(an/an1), in, inf);
    R = sqrt(R);  % 3^(1/2)
    
    d_an = (2*in^2)/(3^in);
    d_an1 = subs(d_an, in, in+1);
    Rd = limit( abs(d_an/d_an1), in, inf);
    Rd = sqrt(Rd); %  3^(1/2)
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    an = (5^in)/( in*(in+1) );
    an1 = subs(an, in, in+1);
    R = limit( abs(an/an1), in, inf);  %1/5
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    an = ( (-2)^in ) / ( in*(in+1)*(in+2) );
    an1 = subs(an, in, in+1);
    R = limit( abs(an/an1), in, inf);
    R = sqrt(R)
end
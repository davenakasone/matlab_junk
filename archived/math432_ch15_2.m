%{
    ch 15.2 power series
    
    #'a'        geo series
    #'b'        exp series
    #6          center and radius of convergence
    #7          " "
    #8          " "
    #9          " "  but don't forget (z-z0)^2n  so R = sqrt(lim)
    #10         " "
    #11         " "  (z-z0)^n term a little fucky
    #12         " "
    #13         " "  watch ^4 on var
    #14         " " ^2
    #15         " " just a lot of factorials
    #17         be careful with var  ^2n+1   that term has to come out
    #18         has a const oo anyway

%}
clc;
close all;
clearvars;


                sel = 9;  % CHANGE CHANGE CHANGE



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
global i_n; syms i_n; assume(i_n, {'real', 'integer'}); % index n
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global N; syms N; assume(N, {'real','integer'});  % bound of series, usually sub with inf

%------------------------------------------------------------------------------------------ #'a'
if sel == 'a'
    z = .5;  % convg if abs < 1 , diverge otherwise...geometric series
    sers = symsum(z^i_n, i_n, 0, N);   % [ 0, N]
    limit(sers, N, inf);  
end


%------------------------------------------------------------------------------------------ #'b'
if sel == 'b'
    z = 1 + 1j;
    seq = z^i_n / factorial(i_n);
    sers = symsum( seq , i_n, 0, N);    % convg for all
    convg = double(limit(sers, N, inf))
    actual = exp(z)
    
    % ratio test
    n1term = subs(seq, i_n, i_n+1);
    ratio = n1term / seq;
    result = limit(ratio, i_n, inf); % should be 0 for all by ratio test  
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    seq = (4^i_n)*(Zt+1)^i_n;      % implies a_n = 4^n  , z0 = -1
    an = 4^i_n;
    an1 = 4^(i_n+1);
    R = limit( abs(an/an1), i_n, inf);
    fprintf(' converges around a circle (-1, 0) r = %.3f\n', R);
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    % z0 implied:  pi/2
    an = ((-1)^i_n)/factorial(2*i_n);
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf) % converges for any z ... R  = oo
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % implies center pi j
    an = i_n^i_n / factorial(i_n);
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf)
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % center j  ...but (z-z0) term is sqared
    an = ( i_n * ( i_n - 1 ) ) / 3^i_n;
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf);
    fprintf('R must be %s because var term is squared\n', sqrt(R));
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    % center 2j  
    an = i_n^i_n;
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf)
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    % z0 implies 0 center
    coeff = (2-1j)/(1+5j); % nothing happening on coeff seq...limit is 1   maybe book ment ^n
    an = coeff^i_n;
    an1 = coeff^(i_n+1);
    R = limit( abs(an/an1), i_n, inf)  % works
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    % 0 center implied
    an = (((-1)^i_n)*i_n)/8^i_n;
    an1 = subs(an,i_n,i_n+1);
    R = limit( abs(an/an1), i_n, inf)
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    % the center implies -j  but it is wrapped up by ^4  so R^(1/4) on the way out
    an = 16^i_n;
    an1 = 16^(i_n+1);
    R = limit( abs(an/an1), i_n, inf);
    rad = R^(1/4)
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    % center 0, sqrt out
    an = ((-1)^(i_n)) / ( (2^(2*i_n)) * (factorial(i_n)^2) );
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf) % R = oo   convg for all Z
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    an = factorial(2*i_n) / ( (4^(i_n)) * factorial(i_n)^2 );
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf)
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    % center 0, but handle the ^2n+1   the +1 is of no consequence
    an = (2^i_n) / (i_n*(i_n+1));
    an1 = subs(an, i_n, i_n+1);
    R = limit( abs(an/an1), i_n, inf);
    rad = sqrt(R) 
end



    
    
    
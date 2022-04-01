%{

    ch 13.6
        
        #1   ex2a
        #2   ex2b
        #3   ex2c
        #4   p4    RC test
        #5   p5    harmonic test
        #6   p6    num eval
        #7   p7    num eval
        #8   p8    num eval
        #9   p9    num eval
        #10  p10   num eval
        #11  p11   find boundary

%}
clc;
close all;
clearvars;

sel = 11;

% Math432 globals   only X и Y are needed in helper funcitons 
global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Z; Z = (X + 1j*Y);

%------------------------------------------------------------------------------------------ #1
if sel == 1
    % cos(z) = 5
    % this means   5 = (1/2) [ exp(jz) + exp(-jz) ]   because cos(z) = (1/2) [exp(jz) + exp(-jz)]
    % this means exp(2jz) - 10exp(jz) + 1 = 0   after multiplying by exp(jz)
    % an quadratic in exp(jz) instead of x
    % by quad formula, exp(jz) = [ 10 +/- sqrt(100 - 4) ] / 2
    % exp(jz) = 5 +/- sqrt(25 - 1)   =>    exp(jz) = 9.899 or .101
    % exp(jz) = exp(-y+jx)   
    % so exp(jx) = 1 and exp(y) is 9.899 or .101
    % but they can solve the eqn in any 2n(pi) cycle
    % so, y = +/- 2.292 , x = 2n(pi)   n = 0, 1, 2,...
    left = cos(Z);
    right = 5;
    %answ = solve(left == right, X, Y); no good
    leftT = (1/2)*( exp(1j*Z) + exp(-1j*Z));
    rightT = 10; % 2 * 5  take (1/2) over
    quad = exp(1j*2*Z) - 10*exp(1j*Z) + 1;
    pol = [1, -10, 1]; % how to solve for roots of any polynomial
    rots = roots(pol);
    rot1 = rots(1);
    rot2 = rots(2);
    
    syms y;
    sol1 = solve( exp(-y) == rot1, y);
    sol2 = solve( exp(-y) == rot2, y);
    display(double(sol2));
    % they are y = +/- 2.2924
    % and exp(jx) has to = 1, so x = 2n(pi)
    %check
    fprintf(' f(z) = cos(z) = cos( x + jy) = 5   ...check\n');
    for k = 1:5
        y =  1j*sol2;
        x = k * 2 * pi;
        z = x + y;
        temp = subs(left, Z, z);
        fprintf('answer = %.4f\n', real(temp))
    end
end


%------------------------------------------------------------------------------------------ #2
if sel == 2  
    % cos(z) = 0
    % turn it into cos(z) = cos(x) cosh(y) - j sin(x) sinh(y)   so x = pi/2 or y = pi
    % that means z =  (pi/2) + n(pi)    n = 0, 1, ...+/-
    % cos(x) = 0 at any  x =pi/2
    % sinh(y) = 0 at y = 0 only 
    % any x = pi/2   (2n + 1) multiple...y stays at 0
    f = cos(Z);
    for k = 1:5
        temp = subs(f, Z, (pi/2)*(2*k+1));
        fprintf(' answer = %s\n',temp);
    end 
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % sin(z) = 0 .....must have sinh(y) = 0, any pi multiple of x to supress sin(x)
    f = sin(Z);
    for k = 1:5
        temp = subs(f, Z, pi*k);
        fprintf('answer = %s\n', temp);
    end
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    % prove cos(z) is entire
    funz1 = cos(Z);
    dump1 = fun_rCRtest(funz1);
    pause(1);
    clc;
    % prove sin(z) is entire
    funz2 = sin(Z);
    dump2 = fun_rCRtest(funz2);
    pause(1);
    clc;
    % prove cosh(z) is entire
    funz3 = cosh(Z);
    dump3 = fun_rCRtest(funz3);
    pause(1);
    clc;
    % prove sinh(z) is entire
    funz4 = cosh(Z);
    dump4 = fun_rCRtest(funz4); 
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    funz1 = imag(cos(Z));
    dump1 = fun_harminic2test(funz1);
    funz2 = real(sin(Z));
    dump2 = fun_harminic2test(funz1);
    %fprintf('%s',funz1);
     
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % find value of sin(j2pi)
    % use  sin(jz) = j sinh(z) = -j sin(jz)
    
    num = 1j*2*sym(pi);
    funz = sin(Z);
    display(num);
    sol = subs(funz, Z, num);
    display(sol);
    sol = double(sol);
    display(sol);
    
    % also
    check = sinh(2*pi);
    check = 1j*check;
    display(check);
end

%------------------------------------------------------------------------------------------ #7
if sel == 7
    % find  cos(j)   ....  cos(jz) = cosh(z)
    num = 1*j;
    funz = cos(Z);
    sol = subs(funz, Z, num);
    display(sol);
    sol = double(sol);
    %display(sol);
    display(cosh(1));  % it is cosh(1)
    
    %find sin(j) ...  sin(jz) = j sinh(z)
    display(sin(1j));
    display(1j*sinh(1));
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % cos(j pi)     since cos(jz) = cosh(z)   ,,,,  cosh(pi) should be equiv
    display(cos(1j*pi));
    display(cosh(pi));
    % cosh(j pi)    since cosh(jz) = cos(z) ,,, cos(pi) should be equiv
    display(cosh(1j*pi));
    display(cos(pi));
end

%------------------------------------------------------------------------------------------ #9
if sel == 9
    % cosh(-1 +j2)    cosh⁡(z1 + z2 )= cosh⁡(z1) cosh⁡(z2) + sinh⁡(z1) sinh⁡(z2)
    % = cosh(-1) cosh(j2) + sinh(-1) sinh(j2)
    % = cosh(-1) cos(2)   + sinh(-1) j sin(2)        cosh(jz) = cos(z), sinh(jz) = j sin(z)
    display(cosh(-1+1j*2));
    z1 = -1;
    z2 = 1j*2;
    equiv = cosh(z1)*cosh(z2) + sinh(z1) * sinh(z2);
    display(equiv);
    equiv = cosh(-1) * cos(2) + sinh(-1) * 1j * sin(2);
    display(equiv);
    pause(1);
    clc;
    
    % cosh(3 + 4j) ... same procedure
    % = cosh(3) cosh(4j) + sinh(3) sinh(4j)
    % = cosh(3) cos(4) + j sinh(3) sin(4)
    display(cosh(3+1j*4));
    display(cosh(3)*cos(4)+1j*sinh(3)*sin(4)); 
end

%------------------------------------------------------------------------------------------ #10
if sel == 10
    % sinh(3 + 4j)      sinh⁡(z1 + z2) = sinh⁡(z1) cosh⁡(z2) + cosh⁡ (z1)  sinh⁡(z2) 
    % = sinh(3) cosh(4j) + cosh(3) sinh(4j)
    % ...sub   =   sinh(3) cos(4) + cosh(3) j sin(4)
    display(sinh(3+1j*4));
    display(sinh(3)*cos(4)+1j*cosh(3)*sin(4));
    display('match');
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    % want sin(z) = 100
    % if sin(z) = ( 1/ 2j ) (  exp(jz)  - exp(-jz) )
    %  100 = ( 1 / 2j ) [ exp( -y + jx)  -  exp ( y + jx) ]
    %  200 j =  exp(-y) exp(jx) - exp(y) exp(jx)    multiply
    % ...or retain z from begninning   200j exp(jz) = exp(2jz) - 1
    pol = [1, -1j*200, -1];  % exp(2jz) - 200 j exp(jz) - 1 = 0
    rotz = roots(pol);
    display(rotz);
    syms z;
    sol1 = solve(exp(z) == rotz(1), z);
    display(double(sol1));
    sol2 = solve(exp(z) == rotz(2), z);
    display(double(sol2));
    % don't forget x is not fixed...can be anything 2pi interval
end

    

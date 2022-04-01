%{
    ch13.4  Cauchy-Riemann criteria, harmonics, 
    
    #2 is    f(z) = jz conj(z)     analytic?
    #3 is    f(z) = exp(-2x) (cos2y - jsin2y)     analytic?
    #4 is    f(z) = exp(x) (cosy - jsiny)     analytic?
    #5 is    f(z) = real(z^2) - j imag(z^2)     analytic? 
    #6 is    f(z) = 1 / (z - z^5)     analytic?
    #7 is    f(z) = j / (z^8)     analytic?            POLAR
    #8 is    f(z) = angle(z pi 2)     analytic? 
    #9 is    f(z) = 3(pi^2) / (z^3 + z4pi^2 )    analytic?        be careful with poles
    #10 is    f(z) = ln(abs(z)) + j angle(z)   analytic?
    #11 is    f(z) = cosx coshy - j sinx sinhy   analytic?
    #12  u = x^2 + y^2      is it harmonic?   find conj = v(x,y) if so
    #13  u = xy      is it harmonic?   find conj = v(x,y) if so     
            fyi    f_int = int(f,x,'IgnoreAnalyticConstraints', true);    so it can provide reasonable result
    #15  u = x / (x^2 + y^2)      is it harmonic?   find conj = v(x,y) if so
    #16  u = sinx coshy      is it harmonic?
%}

clc; 
clf;
close all;
clearvars;
select = 16; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice

%----------------------------------------------------------------------------------------------- 2
if select == 2;     % #2 is    f(z) = jz conj(z)     analytic?  no...  only at (0,0)
    syms x real;
    syms y real;
    z = (x + j*y);
    f = j*z*conj(z);
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 3
if select == 3;     % #3 is    f(z) = exp(-2x) (cos2y - jsin2y)     analytic?   yes, looks good
    syms x real;
    syms y real;
    f = (exp(-2*x)) * ( (cos(2*y)) - j*(sin(2*y)) );
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 4
if select == 4;     % #4 is    f(z) = exp(x) (cosy - jsiny)     analytic?     no
    syms x real;
    syms y real;
    f = (exp(x)) * ( (cos(y)) - j*(sin(y)) );
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 5
if select == 5;     % #5 is    f(z) = real(z^2) - j imag(z^2)     analytic?    no
    syms x real;
    syms y real;
    z = (x + j*y);
    f = (real(z^2)) - j*(imag(z^2));
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 6
if select == 6;     % #6 is    f(z) = 1 / (z - z^5)     analytic?     no, good polar one, one at poles
    syms x real;                         % z(1-z^4) = z (1-z) (1+z) (z-j)( z + j)  are only analytic pts
    syms y real;
    z = (x + j*y);
    f = 1 / (z - (z^5));
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 7
if select == 7;     % #7 is    f(z) = j / (z^8)     analytic?    no,  good polar
    syms r real;
    syms th real;
    x = r*cos(th);
    y = r*sin(th);
    z = (x + j*y);
    f = ( (j*1) / (z^8) );
    u = real(f);
    v = imag(f);
    
    ur = diff(u,r,1);
    uth = diff(u,th,1);
    vr = diff(v,r,1);
    vth = diff(v,th,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ur -(1/r)*vth;
    test2 = vr + (1/r)*uth;
    fprintf("ur = %s     vth = %s     \n", ur, vth);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuth = %s     vr = %s \n', uth, vr);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 8
if select == 8;     % #8 is    f(z) = angle(z pi 2)     analytic?    
    syms x real;                        
    syms y real;
    z = (x + j*y);
    f = angle(2*pi*z);
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 9
if select == 9;     % #9 is    f(z) = 3(pi^2) / (z^3 + z4pi^2 )    analytic?  rationals are good except poles  
    syms x real;                        
    syms y real;
    p = sym(pi); % important to keep accuracy
    z = (x + j*y);
    f = (3*(p^2)) / ( (z^3) + 4*(p^2)*z );
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);    
end


%----------------------------------------------------------------------------------------------- 10
if select == 10;     % #10 is    f(z) = ln(abs(z)) + j angle(z)   analytic?     it is good   z not 0
    % matlab uses log() for natural   and log10() for base 10
    % for negative and complex numbers, use  u = a + j*b  ->  log(abs(u)) + j*angle(u)
    syms x real;                        
    syms y real;
    z = (x + j*y);
    f = log(abs(z)) + j*angle(z);
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 11
if select == 11;     % #11 is    f(z) = cosx coshy - j sinx sinhy   analytic?     simplfy...cosh(y-jx)  good
    syms x real;                        
    syms y real;
    z = (x + j*y);
    f = cos(x) * cosh(y) - j*sin(x)*sinh(y);
    u = real(f);
    v = imag(f);
    
    ux = diff(u,x,1);
    uy = diff(u,y,1);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    
    fprintf("\nf(z) = %s\n\n", expand(f));
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 12
if select == 12     % #12  u = x^2 + y^2      is it harmonic?   find conj = v(x,y) if so
    syms x real;
    syms y real;
    u = x^2 + y^2;
    uxx = diff(u, x, 2);
    uyy = diff(u, y, 2);
    Unamb = uxx + uyy;
    fprintf('\nuxx + uyy = %s   +  %s    =  %s\n', uxx, uyy, Unamb);
    fprintf("it does not solve the Laplace equation...can't be harmonic\n");
end


%----------------------------------------------------------------------------------------------- 13
if select == 13     % #13  u = xy      is it harmonic?   find conj = v(x,y) if so
    syms x real;
    syms y real;
    u = x*y;
    ux = diff(u, x, 1);
    uxx = diff(u, x, 2);
    uy = diff(u, y, 1);
    uyy = diff(u, y, 2);
    Unamb = uxx + uyy;
    fprintf('\nuxx + uyy = %s   +  %s    =  %s\n', uxx, uyy, Unamb);
    fprintf("it solves the Laplace equation...find v(x,y)\n");
    % if ux = vy, integrate ux w/r x   = u + h(y)
    
    v = (1/2)*(y^2)-(1/2)*(x^2);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    f = u + j*v;
    fprintf("\nf(z) = %s\n\n", f);
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2);
end


%----------------------------------------------------------------------------------------------- 15
if select == 15     % #15  u = x / (x^2 + y^2)      is it harmonic?   be careful homo
    syms x real;
    syms y real;
    u = x / (x^2 + y^2);
    ux = diff(u, x, 1);
    uxx = diff(u, x, 2);
    uy = diff(u, y, 1);
    uyy = diff(u, y, 2);
    Unamb = uxx + uyy;
    fprintf('\nuxx + uyy = %s   +  %s    =  %s\n', uxx, uyy, Unamb);
    fprintf("it actually works, just need some algebra\n");
    
    v = -int(uy,x);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    f = u + j*v;
    fprintf("\nf(z) = %s\n\n", f);
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2)
end


%----------------------------------------------------------------------------------------------- 16
if select == 16     % #16  u = sinx coshy      is it harmonic?   
    syms x real;
    syms y real;
    u = sin(x)*cosh(y);
    ux = diff(u, x, 1);
    uxx = diff(u, x, 2);
    uy = diff(u, y, 1);
    uyy = diff(u, y, 2);
    Unamb = uxx + uyy;
    fprintf('\nuxx + uyy = %s   +  %s    =  %s\n', uxx, uyy, Unamb);
    fprintf("it actually works, just need some algebra\n");
    
    v = cos(x) * sinh(y);
    vx = diff(v,x,1);
    vy = diff(v,y,1);
    f = u + j*v;
    fprintf("\nf(z) = %s\n\n", f);
    test1 = ux - vy;
    test2 = uy + vx;
    fprintf("ux = %s     vy = %s     \n", ux, vy);
    fprintf("test1 says: %s\n",test1);
    fprintf('\nuy = %s     vx = %s \n', uy, vx);
    fprintf('test2 says: %s \n', test2)
end

%{
    pro tip:
                z = -1 + 1j;
    display(rad2deg(atan2( imag(z), real(z))));    % keep quadrant so angle is always good


    ch13.5 , exponential function f(z) = exp(z)
    
    #1  ex1
    #2  ex2
    #3  p2
    #4  p3
    #5  p4
    #6  p5
    #7  p6
    #8  p7
    #9  p8  exp form
    #10 p9
    #11 p10  sqrt(j) , sqrt(-j)
    #12 p11  
    #13 p12   1/(1-z)   wtf?
    #14 p13
    #15 p14
    #16 p15
    #17 p16
    #18 p17
    #19 p18a   using fun_rCRtest()  looks good
    #20 p18b                                        GRAPHs...contour и surface

%}
%color = uisetcolor([1, 1, 0], 'Selecf Color'); % .9, .9, .9 is nice
%clear all;  if things aren't going away
%consts = cls_CONST();  % common physical constants
%consts.help();
%ee = cls_EE330_helper(); % vectors, cyn, cph, rec, 

clc;
close all;
clearvars;


sel = 20;

% try this out
global X; syms X; assume(X, 'real');
global Y; syms Y; assume(Y, 'real');
global Z; Z = (X + 1j*Y);

%------------------------------------------------------------------------------------------ #1
if sel == 1
    syms z;
    f = exp(z);
    x = 1.4;
    y = -.6;
    val = (1.4 - 1j*.6);
    result = subs(f, z, val);
    display(double(result));
    % note:
    eqiv = exp(x+1j*y);
    display(eqiv);
    mag = abs(eqiv);
    display(mag);
    angl = angle(eqiv);
    display(angl);
end

%------------------------------------------------------------------------------------------ #2
if sel == 2
    syms x real;
    syms y real;
    z = x + 1j*y;
    f = exp(z);
    target = 3 + 1j*4;
    %answ = solve(f == target, x, y);
    %display(target);                       solve f(z) = exp(z) = 3 + j4
    mag = abs(target);
    fprintf('the abs is %.4f\n',mag);
    a = log(5);
    fprintf('abs = 5  means abs(exp(z)) = 5 = exp(x)  .... so x = ln(5) : %.4f\n', a);
    b = acos(3/mag);
    fprintf('focus on real part of Euler eqn, exp(z) = exp(x)[ cosy + jsiny]\n');
    fprintf(' 3 = exp(x) cos(y)  for real part... y = acos(3 / ln(5) ) = %.4f \n', b);
    % or exp(x) siny = 4   either way
    check = double(subs(f, z, a+1j*b));
    display(check);  % actually inf sol...fixed x, y varies 2pi n
end

%------------------------------------------------------------------------------------------ #3
if sel == 3
    target = 3 + 1j*4;
    syms z;
    f = exp(z);
    mag = abs(target);
    y = acos( real(target) / mag);
    x = log(mag);
    check = x + 1j*y;
    sol = subs(f, z, check);
    display(double(sol));
end

%------------------------------------------------------------------------------------------ #4
if sel == 4
    target = 2*1j*pi*(1+1*j);
    a = real(target);
    b = imag(target);
    mag = abs(target);
    x = log(mag);
    y = acos( a / mag );
    check = x + 1j*y;
    syms z;
    f = exp(z);
    sol = double(subs(f, z, check));
    display(double(target));
    display(sol);
end

%------------------------------------------------------------------------------------------ #5
if sel == 5
    target = .6 - 1j*1.8;
    mag = abs(target);
    x = log(mag);
    y = asin( imag(target) / mag);
    check = x + 1j*y;
    syms z;
    f = exp(z);
    display(target);
    display(double(subs(f, z, check)));
end

%------------------------------------------------------------------------------------------ #6
if sel == 6
    tgt = 2 + 1j*pi*3;
    mag = abs(tgt);
    x = log(mag);
    y = acos(2/mag);
    chk = x + 1j*y;
    syms z;
    f = exp(z);
    disp(tgt);
    disp(double(subs(f, z, chk)));
end

%------------------------------------------------------------------------------------------ #7
if sel == 7
    tgt = 1j*11*pi / 2;     % the target is completley imaginary  --> y = 90°
    mag = abs(tgt);
    x = log(mag);
    %y = asin( imag(tgt) / mag); % have to use sin...no imag part..or just see angle
    y = acos( real(tgt) / mag);
    check = x + 1j*y;
    syms z;
    f = exp(z);
    answ = double(subs(f, z, check));
    disp(tgt);
    disp(answ);      % notice the negative zero convention....not same as IEEE convention
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    tgt = sqrt(2) + 1j*(1/2)*pi;
    mag = abs(tgt);
    x = log(mag);
    y = acos( real(tgt) / mag);
    chk = x + 1j*y;
    syms z;
    f = exp(z);
    sol = double(subs(f,z,chk));
    display(double(tgt));
    display(sol);
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % want exp(z^(1/n)) in exp form
    %  z^(1/n) = { r [ cos(t) + j sin(t) ] } ^(1/n)
    %  z^(1/n) = { r^(1/n) [ cos(2kpi + t / n) + j sin(cos(2kpi + t / n) ]
    % so exp(z) = r^(1/n) * exp( j (2kpi + t / n))    k = 0, 1, ..., n-1
    syms n real;
    syms k real;
    syms r real;
    syms thta real;
    syms z;
    f = exp(z);
    display(f);
    mag = r^(1/n);
    th = 2*k*pi + thta / n;
    g = mag * subs(f,z,th);
    display(g);
end

%------------------------------------------------------------------------------------------ #10
if sel == 10
    tgt = 4 + 1j*3;
    th = angle(tgt);
    mag = abs(tgt);
    f = mag * exp(1j*th);
    fprintf(' f(z) = %.1f * exp(j * atan(%.1f / %.1f) ) \n', mag, imag(tgt), real(tgt));
    fprintf(' f(z) = %.1f * exp(j * %.2f°) \n', mag, rad2deg(th));
    display(f);
end

%------------------------------------------------------------------------------------------ #11
if sel == 11
    tgt = 0 + sqrt(1j);
    mag = abs(tgt);
    th = pi/2;  % implied
    % Euler > z = r [ cos(t) + j sin(t) ] > r exp(jt)
    % look at sqrt(j) in polar >  1 * [ cos(pi/2) + j sin(pi/2) ] ^ (1/2)
    % demovire > = cos(2pk + p/2 / 2) + j sin( 2pk + p/2 / 2 )   k = 0, 1
        % = cos( (4k + 1)p / 4 ) + j sin( (4k + 1)p / 4 )   k = 0, 1
        
    % similar for sqrt(-j)   angle is -pi/2   abs = 1
    % z = sqrt(-j) > 1 * [ cos(-pi/2) + j sin(-pi/2) ] ^ (1/2)
    % >> [ cos(pi/2) - j sin(pi/2) ] ^ (1/2)         cos even, sin odd
    %   demovire
    
        fprintf('sqrt(j) = exp( (4k + 1)pi / 4 )    k = 0, 1\n');
        fprintf('sqrt(-j) = exp( -(4k + 1)pi / 4 )  k = 0, 1\n');
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    tgt = -6.3;  % the magnitude is real, angle is 180°
    % f = 6.3 * exp(j*pi)
    display(6.3*exp(1j*pi));
end

%------------------------------------------------------------------------------------------ #13
if sel == 13
    % you want to write f(z) = 1 / (1-z) in exp form
    % the polar form of z is   z = r ( cos(t) + jsin(t) )
    % if z = x + jy in rec, r = sqrt(x^2 + y^2) , t = atan(y/x)   ! t = atan2( imag(z), real(z) )
    % 1 / ( 1-z) = 1 / ( 1 - r exp(jt) )
    % multiply to and bottom -> trig sub
    %   use multiplication by   (1 - r exp(-jt)) / (1 - r exp(-jt))
    fprintf(' f(z) = [ 1 - r exp(-jt) ] / [ 1 + r^2 - 2 r cos(t) ] in exponential\n');
end

%------------------------------------------------------------------------------------------ #14
if sel == 14
    tgt = 1 + 1j;
    th = atan2( imag(tgt), real(tgt) );
    r = abs(tgt);
    fprintf(' f(%.1f + j%.1f) = %.4f * exp( j * %.1f°) \n',...
        real(tgt), imag(tgt), r, rad2deg(th));
end

%------------------------------------------------------------------------------------------ #15
if sel == 15
    % find real and imaginary parts of exp(-pi z)
    % assume z = x + jy
    % -pi z = -pi (x + jy)      so exp(-pi z) = exp( -pi x) exp(-pi j y)
    % >> exp( -pi x ) [ cos(-pi y) + j sin(-pi y) ]
        %   exp( -pi x ) [ cos(pi y) - j sin(pi y) ]
            % real =>  exp(-pi x) cos( pi y)
            % imag =>  - exp(-pi x) sin( pi y)
            display('see comments in #15');
    syms x real;
    syms y real;
    z = x + 1j*y;
    f = exp(-pi*z);
    display(real(f));
    display(imag(f));
end
 
%------------------------------------------------------------------------------------------ #16
if sel == 16
    % exp(z^2)
    % if z = x + jy    --> z^2 = x^2 - y ^2 + j 2 x y
    % f( exp(z^2) = [exp(x^2) / exp(y^2) ] * exp(j2xy)
        %  f  =  [exp(x^2) / exp(y^2) ]  *  { cos(2xy) + j sin(sxy) }
            % real = [exp(x^2) / exp(y^2) ] * cos(2xy)
            % imag = [exp(x^2) / exp(y^2) ] * sin(2xy)
            display('see comments in #16');
     syms x real;
     syms y real;
     z = x + 1j*y;
     f = exp(z^2);
     display(real(f));
     display(imag(f));
end

%------------------------------------------------------------------------------------------ #17
if sel == 17
    % f( exp(1/z) )
    % z = x + jy       --> 1 / z = 1 / x + jy
    % multiply by conj   1/z --> x - jy / x^2 + y^2   --> x / x^2 + y^2  - jy / x^2 + y^2
    % 
    display('see comments in #17');
    syms x real;
    syms y real;
    z = x + 1j*y;
    f1 = 1/z * (1/conj(z));
    display(simplify(f1));
    f2 = conj(z) * simplify(f1);
    display(f2);
    g = exp(f2);
    gr = real(g);
    gi = imag(g);
    display(gr);
    display(gi);
end

%------------------------------------------------------------------------------------------ #18
if sel == 18
    % want exp(z^3)
    syms x real;
    syms y real;
    z = x + 1j*y;
    f = exp(z^3);
    display(simplify(real(f)));
    display(simplify(imag(f)));
end

%------------------------------------------------------------------------------------------ #19
if sel == 19
    f = Z^2; % testing this out
    bol = fun_rCRtest(f); % it looks good
    if bol == 1
        fprintf('f(z) = %s passed CR test\n', f);
    else
        fprintf('f(z) = %s failed\n', f);
    end
    
    fprintf('\n\n');
    f = exp(Z); 
    bol = fun_rCRtest(f); 
    
    fprintf('\n\n');
    f = exp(1/Z); 
    bol = fun_rCRtest(f);
    
    fprintf('\n\n');
    f = exp(conj(Z)); 
    bol = fun_rCRtest(f); 
    
    fprintf('\n\n');
    syms k;
    f = exp(X) * (cos(k*Y) + 1j*sin(k*Y));
    bol = fun_rCRtest(f);  
end

%------------------------------------------------------------------------------------------ #20
if sel == 20
    % f(z) = exp(z) = exp(x + jy) is real when y is 0  exp(x) [ cos(y) + j sin(y) ]
    %fun = (Z)^2;
    %toss = fun_graphSurCon(fun);
    
    %fun = exp(Z);
    %toss = fun_graphSurCon(fun);
    
    %num = 3 - 1j*4;
    %rots = fun_demovNum(num, 1, 4);
    %display(rots);
    
    %tgt = exp(1/Z);
    %fun_graphSurCon(tgt);
    
    %fun = exp(Z);
    %fun_tetraView(fun, -5, -5);
    
    %fun = exp(-Z);
    %fun_graphSurCon(fun);
    
    %fun1 = exp(conj(Z));
    %fun2 = conj(exp(Z));
    %fun_graphSurCon(fun1);
    %fun_tetraView(fun1, 0, 0)
    %fun_graphSurCon2(fun1,fun2);
    
    funz = Z^2;
    %fun_tetraView(funz, 0, 0);
    fun_graphSurCon(funz);
end
    

    
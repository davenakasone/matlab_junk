%{
    chapter 2 : frequency domain, using Laplace

    1  :  appxB1, basic complex >> real(), imag(), angle(), atan2(), abs()
    2  :  appxB2, represent a polynomial by coefs
    3  :  appxB3, poly() vs root() vs coef()
    4  :  appxB4, represent polynomial, roots known >> poly(), syms, expand(), pretty(), size(), for...
    5  :  appxB5, find the roots of a given polynomial >> roots()
    6  :  appxB6, poly multiplication by conv()
    7  :  appxB7, get PFE  ... even does poly division if needed >> residue(), sym2poly()
    8  :  appxB8, PFE ex, use that fliplr() to get right...
    9  :  appxB9, make transfer function...just use >> tf(), zpk()
    10 :  appxB10, convert TF from poly to zero/pole K >> tf2zp(), zp2tf()
    11 :  try g_ct_2d()
    12 :  appxF1, Laplace and inverse >> laplace(), ilaplace()
    13 :  appxF2, >> collect(), expand(), factor(), simple(), simplify(), vpa()
    14 :  appxF3, isolate a TF >> numden() 
    15 :  appxF4, use the matrix math >> det(),  \
    16 :  TryIt2.2, basic laplace()
    17 :  TryIt2.3/4, basic tf()
    18 :  ex2.2, laplace using PFE, rat() to make fraction
    19 :  ex2.5, from transfer function
    20 : hw 2.41, linearize
  

my_coeffs = fliplr(coeffs(my_poly, var, 'All'));
be able to manipulate common functions
be able to see a 2D graph
%}
close all;
clc;
select = 20;


%------------------------------------------------------------------------------------------
if (select == 1)
    mm = 2 + 2j;
    fprintf("\nthis is a complex number:  %d\n", mm)
    mm_rad = angle(mm);
    mm_deg = rad2deg(mm_rad);
    mm_mag = abs(mm);
    fprintf("\tphase:  %0.2f [rad]  %0.2f [deg]  ,  mag:  %0.2f\n",...
        mm_rad, mm_deg, mm_mag);
    % another way to get angle
    mm_ang = atan2(imag(mm), real(mm));
    fprintf("\tdon't for get about atan2()  ,  radians:  %0.2f\n", mm_ang);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    my_poly = [1, 7, -3, 25] 
    fprintf("\nthis is: {%d} x^3 + {%d} x^2 + {%d} x^1 + {%d} x^0\n",...
        my_poly(1), my_poly(2), my_poly(3), my_poly(4));
    fprintf("\tremeber matlab indexs on (), not [], and they start at 1, not 0\n");
end


%------------------------------------------------------------------------------------------
if (select == 3)
    my_poly = [-1, -1]    % (x+1)^2  ... you need the roots
    my_roots = roots(my_poly) % -1, -1 as roots , confirmed
    
    syms x;
    my_expp = (x - my_roots) .* (x - my_roots);
    pretty(my_expp)
    pretty(expand(my_expp))
    my_coeffs = coeffs(my_expp)
end


%------------------------------------------------------------------------------------------
if (select == 4)
    my_roots = [-2, -5, -6];    % (x+2)(x+5)(x+6)
    store = 1;
    syms x;
    
    sizz = size(my_roots);    % elements in (rows, cols)
    sizzr = sizz(1);          % number of elements in row
    sizzc = sizz(2);          % number of elements in col
    
    for ii = 1:1:sizzc
        store = store .* (x - my_roots(ii));
    end
    
    fprintf("\nresulting polynomial...factored and expanded\n");
    pretty(store);
    pretty(expand(store)); 
end


%------------------------------------------------------------------------------------------
if (select == 5)
    my_poly = [5, 7, 9, -3, 2];    % x^4  ... x^0 , just coeffiecients needed
    my_roots = roots(my_poly)
end


%------------------------------------------------------------------------------------------
if (select == 6)
    poly1 = [1, 7, 10, 9]
    poly2 = [1, -3, 6, 2, 1]
    poly12 = conv(poly1, poly2);
    
    syms x;
    collect = 1;
    for ii = 1:1:size(poly12, 2)
        collect = collect .* (x - poly12(1, ii));
    end
    fprintf("\nusing conv() to multiply 2 polynomials:\n");
    pretty(expand(collect));
end


%------------------------------------------------------------------------------------------
if (select == 7) % look at the combination of a vector and poly...
    syms s;
    Fsn = 7*s^2 + 9*s + 12;
    Fsd = s * (s+7) * (s^2 + 10*s + 100);
    Fs = Fsn / Fsd;
    numF = [7, 9, 12];
    denF = conv(poly([0, -7]), [1, 10, 100]);
    [rr, pp, kk] = residue(numF, denF)
    
    fprintf("{r} is the residue that goes in PFE numerator\n");
    fprintf("{p} is root of denominator, all denominators are degree 1\n");
    fprintf("{k} is the constant if there was any poly division\n");
    check = double(subs(Fs, s, pp(1)))% wtf...
    
    [rr, pp, kk] = residue(sym2poly(Fsn), sym2poly(Fsd)) % best...
end


%------------------------------------------------------------------------------------------
if (select == 8)
    syms s;
    numF = 32;
    denF = poly2sym(fliplr([0, -4, 8]), s)
    [rr, pp, kk] = residue(numF, sym2poly(denF));
end


%------------------------------------------------------------------------------------------
if (select == 9) % makes a LTI object
    syms s;
    % TF by vectors in poly form
    Fsn = 150 * (s^2 + 2*s + 7);
    Fsd = s * (s^2 + 5*s + 4);
    Fs = Fsn / Fsd;
    Fsn_c = sym2poly(Fsn)
    Fsd_c = sym2poly(Fsd)
    Fs_t = tf(Fsn_c, Fsd_c)
    
    % TF by vectors in factored form
    Fsn = (s + 2) * (s + 4); % must find zeros
    Fsd = (s + 7) * (s + 8) * (s + 9); % must find poles
    kk = 20; % gain, must be const
    zz = roots(sym2poly(Fsn)); % ... get the roots for poles and zeros
    pp = roots(sym2poly(Fsd));
    Fs_t = zpk(zz, pp, kk);
    
    % TF by LTI object...no
    z = tf('z');
    Fzn = 150 * (z^2 + 2*z + 7);
    Fzd = z * (z^2 + 5*z + 4);
    Fz = Fzn / Fzd
end


%------------------------------------------------------------------------------------------
if (select == 10)
    syms s;
    Fsn = 10*s^2 + 40*s + 60;
    Fsd = s^3 + 4*s^2 + 5*s + 7;
    Fs = Fsn / Fsd;
    Fsn_c = sym2poly(Fsn);
    Fsd_c = sym2poly(Fsd);
    [Fsn_z, Fsd_p] = tf2zp(Fsn_c, Fsd_c)
end


%------------------------------------------------------------------------------------------
if (select == 11)
    syms s;
    start = -6;
    stop = 6;
    dots = 100;
    push_s = linspace(start, stop, dots);
    
    Fs = cos(s);
    push_F = double(subs(Fs, s, push_s));
    Gs = sin(s);
    push_G = double(subs(Gs, s, push_s));
    Hs = Fs * Gs * s^2;
    push_H = double(subs(Hs, s, push_s));
    g_ct_2d(push_s, [push_F; push_G; push_H]); 
end


%------------------------------------------------------------------------------------------
if (select == 12)
    syms s;
    syms t;
    Fs = 2 / ((s+1)*(s+2)^2);
    pretty(Fs);
    Ft = ilaplace(Fs, s, t);
    pretty(Ft);
    Fss = laplace(Ft, t, s);
    pretty(Fss);
    
    indvar = linspace(0.0001, 10, 100);
    gFs = subs(Fs, s, indvar);
    gFt = subs(Ft, t, indvar);
    g_ct_2d(indvar, [gFt; gFs]);
end


%------------------------------------------------------------------------------------------
if (select == 13)
    syms s;
    syms t;
    Ft = 2*exp(-t) - s*t*exp(-2*t) - 2*exp(-2*t);
    pretty(Ft);
    Fs = laplace(Ft, t, s);
    pretty(Fs);
    pretty(simplify(Fs)); % combines
    
    Ft = (3/5) - (1/3)*exp(-4*t);
    Fs = laplace(Ft, t, s);
    pretty(vpa(Fs, 3));
end


%------------------------------------------------------------------------------------------
if (select == 14)
    syms s;
    Fs = 54*(s+27)*(s^3 + 52*s^2 + 37*s + 73) / (s*(s^4 + 872*s^3 + 437*s^2 + 89*s + 65));
    pretty(Fs);
    [numF, denF] = numden(Fs);
    numFp = sym2poly(numF);
    denFp = sym2poly(denF);
    F_tf = tf(numFp, denFp)
    G_zpk = zpk(F_tf)
end


%------------------------------------------------------------------------------------------
if (select == 15)
    syms s;
    syms R1;
    syms R2;
    syms L;
    syms C;
    syms V;
    
    A2 = [R1+L*s, V; -L * s, 0];
    A = [R1+L*s, -L*s; -L*s, 1/(C*s)];
    
    I2 = simplify(det(A2) / det(A));
    pretty(I2);
end


%------------------------------------------------------------------------------------------
if (select == 16)
    syms t;
    syms s;
    Ft = 2*exp(-t) - 2*t*exp(-2*t) - 2*exp(-2*t);
    Fs = laplace(Ft, t, s);
    pretty(Fs);
end


%------------------------------------------------------------------------------------------
if (select == 17)
    syms s;
    Fs = 3/(s*(s^2 + 2*s +5));
    [numF, denF] = numden(Fs);
    F_tf = tf(sym2poly(numF), sym2poly(denF))
    F_t = ilaplace(Fs);
    pretty(F_t);
end


%------------------------------------------------------------------------------------------
if (select == 18)
    syms t;
    syms s;
    F_s = 10 / (s * (s+2) * (s+3)^2);
    
    [numF, denF] = numden(F_s);
    [rr, pp, kk] = residue(sym2poly(numF), sym2poly(denF));
    rat(rr)% shows pfe
    
    f_t = ilaplace(F_s);
    pretty(f_t);
end


%------------------------------------------------------------------------------------------
if (select == 19)
    syms s;
    syms t;
    G_s = s / ((s+4)*(s+8));
    R_s = 1/s^2; % the ramp input
    C_s = G_s * R_s; % because input*tf = output
    c_t = ilaplace(C_s, s, t);
    pretty(c_t)
end


%------------------------------------------------------------------------------------------
if (select == 20)
    syms s;
    syms t;
    Xs_n = 3;
    Xs_d = s * (s^3 + 10*s^2 + 20*s + 30);
    Xs = Xs_n / Xs_d;
    rootz = roots([1, 10, 20, 30])
  
    [rr, pp, kk] = residue(3, [1, 10, 20, 30]);
    x_t = ilaplace(Xs, s, t);
    pretty(x_t)
end
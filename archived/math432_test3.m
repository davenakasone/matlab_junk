%{
    Math_432, test 3

    #1
    #2
    #3
    #4
    #5
    #6
    #7
    #8
    #9
    #10
    #11
    #12
    #13
    #14
    #15
    #16
    #17
    #18
    #19
    #20

    #714 ch17, r14, w=z^2, y[0,2] any x
    #721 ch17, r21, w=1/z, (x-1/2)^2 + y^2 = (1/2)^2  y>0
    #723 ch17, r23, LFT map
    #725 ch17, r25, LFT map
    #726 ch17, r26, LFT map
    #728 ch17, r27, LFT map

    #911 review, 1a residue @ z0 = 0
    #912 review, 1b residue @ z0 = 0
    #913 review, 1c residue @ z0 = 0
    #914 review, 1d residue @ z0 = 0
    #915 review, 1e residue @ z0 = 0
    #916 review, 1f residue @ z0 = 0
    #917 review, 1g residue @ z0 = 0
    #918 review, 2a residue @ pole
    #919 review, 2b residue @ pole
    #920 review, 3a int by residue
    #921 review, 3b int by residue
    #922 review, 3c int by residue
    #923 review, 3d int by residue
    #924 review, 3e int by residue
    #925 review, 3f int by residue
    #926 review, 3g int by residue
    #927 review, 3h int by residue
    #928 review, 3i int by residue
    #929 review, 3j int by residue
    #930 review, 3k int by residue
    #931 review, 3l int by residue
    #932 review, 4a fixed points
    #933 review, 4b fixed points
    #934 review, 5  map
    #935 review, 6  map
    #936 review, 7  map
    #937 review, 8  map
    #938 review, 9  LFT/mobious
    #939 


coeffs(f,var)  makes it lowest degree (left) to highest degree (right)

%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
format rational; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]

global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Zxy; syms Zxy; Zxy = X + 1j*Y;% just use it as a hold and imply  when needed
global Z; syms Z;
global U; syms U; assume(U, 'real');   % Ur ( X, Y)  ...real part of Z = U(X,Y) + j V(X,Y)      
global V; syms V; assume(V, 'real');   % Vi ( X, Y)  ...imag part of Z = U(X,Y) + j V(X,Y)
global W; syms W;

                                select = 7;  % CHANGE CHANGE CHANGE

mat = cls_math('Themistocles');
%------------------------------------------------------------------------------------------ #1
if select == 1
    f = Z^2/(5*Z^3 - 3j*Z^2);
    df = diff(f,Z,1);
    pretty(simplify(df));
    tst = 1/(5*Z-3j);
    d_tst = (-5)*(5*Z-3j)^-2;
    chk = simplify(d_tst-df)
    
end


%------------------------------------------------------------------------------------------ #2
if select == 2
    
end


%------------------------------------------------------------------------------------------ #3
if select == 3
    
end


%------------------------------------------------------------------------------------------ #4
if select == 4
    syms t;
    intg = (t+(4-t^2))*(1-2j*t);
    pretty(expand(intg));
    intgg = int(intg,t);
    pretty(intgg);
    ul = subs(intgg,t,1)
    ll = subs(intgg,t,-2)
    answww = ul-ll
    
end


%------------------------------------------------------------------------------------------ #5
if select == 5
    f = cosh(Z);
    df = diff(f,Z,1);
    ddf = diff(df,Z,1);
    cosh(0)
    
end


%------------------------------------------------------------------------------------------ #6
if select == 6
    
end


%------------------------------------------------------------------------------------------ #7
if select == 7
    syms n;
    an = ((1+3j)^n)/(2^(3*n+1));
    sss = symsum(an, n, 1, inf)
    
    chk = (1+3j)/(2^4);
    chkk = 1/(1-( (1+3j)/(8) ) );
    chkkk = chk*chkk
    
    %{
    rgx = [1,10];
    rgy = [1,10];
    w = ( (Z-1-1j)^2 -1j )/( (Z-1-1j)^2 +1j);
    %mat.fun_image_z(tran, zformat, rng_rx, rng_ty)
    mat.fun_image_z(w, 1, rgx, rgy);
    %}
    
end


%------------------------------------------------------------------------------------------ #8
if select == 8
    rgx = [0,10];
    rgy = [-1,1];
    w = ( sin(pi*1j*Z/2) -1j)/( sin(pi*1j*Z/2) + 1j);
    %mat.fun_image_z(tran, zformat, rng_rx, rng_ty)
    mat.fun_image_z(w, 1, rgx, rgy);
    
end


%------------------------------------------------------------------------------------------ #9
if select == 9
    rgx = [0,pi/6];
    rgy = [0,10];
    w = ( Z^6 -1j)/( Z^6 + 1j);
    %mat.fun_image_z(tran, zformat, rng_rx, rng_ty)
    mat.fun_image_z(w, 2, rgy, rgx);
    
end


%------------------------------------------------------------------------------------------ #10
if select == 10
    
end


%------------------------------------------------------------------------------------------ #11
if select == 11
    
end


%------------------------------------------------------------------------------------------ #12
if select == 12
    
end


%------------------------------------------------------------------------------------------ #13
if select == 13
    
end


%------------------------------------------------------------------------------------------ #14
if select == 14
    
end


%------------------------------------------------------------------------------------------ #15
if select == 15
    
end


%------------------------------------------------------------------------------------------ #16
if select == 16
    
end


%------------------------------------------------------------------------------------------ #17
if select == 17
    
end


%------------------------------------------------------------------------------------------ #18
if select == 18
    
end


%------------------------------------------------------------------------------------------ #19
if select == 19
    
end


%------------------------------------------------------------------------------------------ #20
if select == 20
    
end



%******************************************************************************************
%******************************************************************************************
%******************************************************************************************



%------------------------------------------------------------------------------------------ #911
if select == 911
    f = (Z^2 -1j*2*Z+7)/Z^3;
    pfe = partfrac(f,Z,'FactorMode', 'full');
    % .... just pick that residue off the z^-1 term
end


%------------------------------------------------------------------------------------------ #912
if select == 912
    f = ( (Z+2)*sin(Z) ) / Z^4;
    pfe = partfrac(f,Z,'FactorMode', 'full');
end


%------------------------------------------------------------------------------------------ #915
if select == 915
    f = ( exp(5*Z)-1 ) / (sin(Z))^2 ;
    pfe = partfrac(f,Z,'FactorMode', 'full');
end


%------------------------------------------------------------------------------------------ #919
if select == 919
    f = cos(sym(pi)*1j*Z)/(Z+1j)^2;
    df = diff(f,Z,1);
    pretty(rewrite(simplify(df),'sincos'));
    res = subs(df,Z,1j)
end


%------------------------------------------------------------------------------------------ #714
if select == 714
    w=Z^2;
    mat.fun_image_z(w, 2,[0,1], [0,2*pi]);
end


%------------------------------------------------------------------------------------------ #721
if select == 721
    f = 1/Z;
    w = U + 1j*V;
    wuv = subs(f,Z,w);
    x = real(wuv);
    y = imag(wuv);
    tx = (x-(1/2))^2;
    ty = y^2;
    alph = (U^2 + V^2);
    num = expand((2*U-alph)^2);
    expand(alph^2);
    tt = simplify(num+4*V^2 -alph^2)
end


%------------------------------------------------------------------------------------------ #723
if select == 723
    zz = [-1,0,1];
    ww = [4+3j, 5j/2, 4-3j];
    %{
    tst = mat.fun_lft(zz, ww);
    t1 = (8+1j)*(-4+3j)*(Z+1);
    t2 = (-8+11j)*(4+3j)*(Z-1);
    t3 = t1+t2;
    t4 = (Z-1)*(11j-8)-(Z+1)*(8+1j);
    test1 = subs(tst, Z, zz(1));
    test2 = subs(tst, Z, zz(2));
    test3 = subs(tst, Z, zz(3));
    %}
    
    tw = ( (ww(2)-ww(3))/(ww(2)-ww(1)) );
    tz = ( (Z+1)/(1-Z) ) * (1/tw);
    tzz = (-4+3j)*tz+(4+3j);
    W = tzz/(1-tz);
    subs(W, Z, 1)
end


%------------------------------------------------------------------------------------------ #725
if select == 725
    zz = [1,1j,-1j];
    ww = [1j, -1, 1];
    tst = mat.fun_lft(zz, ww);
    tw = ((ww(2)-ww(3))/(ww(2)-ww(1)));
    tz = ((zz(2)-zz(3))/(zz(2)-zz(1)));
    t1 = subs(tst,Z,zz(1))
    %t2 = subs(tst,Z,zz(2))
    %t2 = subs(tst,Z,zz(3))
end


%------------------------------------------------------------------------------------------ #726
if select == 726
    zz = [0,1,2];
    ww = [2j, 1+2j, 2+2j];
    tst = mat.fun_lft(zz, ww);
    tw = ((ww(2)-ww(3))/(ww(2)-ww(1)));
    tz = ((zz(2)-zz(3))/(zz(2)-zz(1)));
    t1 = subs(tst,Z,zz(1))
    t2 = subs(tst,Z,zz(2))
end


%------------------------------------------------------------------------------------------ #728
if select == 728
    zz = [-1,-1j,1j];
    ww = [1-1j, 2, 0];
    tst = mat.fun_lft(zz, ww);
    tw = ((ww(2)-ww(3))/(ww(2)-ww(1)));
    tz = ((zz(2)-zz(3))/(zz(2)-zz(1)));
    t1 = subs(tst,Z,zz(1))
    t2 = subs(tst,Z,zz(2))
end
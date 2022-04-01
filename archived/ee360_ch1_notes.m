%{
    chapter 1 signals    notes and problems

    #'a'      ex1.1       time shift, reverse, scale  continous    see fun_conT() and fun_conTabg()
    #'b'      discretes   see fun_disN()
    #'c'      unit impulse and step disc
    #'d'      unit impulse and step cont
    #1        pol 2 cart
    #2        cart 2 pol
    #3        limits()  @ inf    be careful with those descrete limits
    #4        compare shifts, reversals, ect   dis
    #5        " " cont
    #6        determine intervals
    #21         HW #1  p21
    #22         HW #2  p22 a-f
    #49         HW #4  a-g, i, k, l
    #55         HW #6  a-e


%}
clc;
close all;
clearvars;


                sel = 55;  % CHANGE CHANGE CHANGE
                
                
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var
global in; syms in; assume(in, {'real', 'integer'}); % index n
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
        
%------------------------------------------------------------------------------------------ #'a'
if sel == 'a'
    imx = [-3,3];
    ptList = [90,90]; %[0,0;1,1;2,1];
    myFun = piecewise(pt<=0, 0, 0<pt<=1, 1, 1<pt<2, 2-pt, pt>=2, 0);
        %fun_conT(myFun, imx, ptList);  % original
    % x(t) at t = t_0    occurs in x(t+1) at  t = t_0 - 1
    % signal advance 1 unit of time
    al = 1;
    be = 1;
    ga = 1;
        %fun_conTabg(myFun, al, be, ga, imx, [900,900]);
    % advance and reverse
    al = -1;
    be = 1;
        %fun_conTabg(myFun, al, be, ga, imx, dmx, ptList);
    % compress horizontal
    al = 3/2;
    be = 0;
        %fun_conTabg(myFun, al, be, ga, imx, dmx, ptList);
    % compress horizontal and advance
    al = 3/2;
    be = 1;
        %fun_conTabg(myFun, al, be, ga, imx, dmx, ptList);
end


%------------------------------------------------------------------------------------------ #'b'
if sel == 'b'
    myFun = i_n^3;
    range = [-3, 5];
    ptl = [0,0; 3,1];
    %fun_disN(myFun, range, [90,90]);
    funA = i_n^2;
    Arng = [3,9];
    funB = i_n^3; %sin(i_n)*exp(i_n/2);
    Brng = [-4,3];
    fun_disN2(funA, Arng, funB, Brng, 0)
end


%------------------------------------------------------------------------------------------ #'c'
if sel == 'c'
    uImp = piecewise(i_n == 0, 1, 0);
    rng = [-5,5];
        %fun_disN(uImp, rng, 2);
            fun_disN(cos(i_n), rng, 2); % stem vs quiv
    
    uStp = piecewise(i_n>=0, 1, 0);
        %fun_disN(uStp, rng, 3);
end


%------------------------------------------------------------------------------------------ #'d'
if sel == 'd'
    uS = piecewise(pt>=0, 1, 0);
    rng = [-5,5];
        %fun_conT(uS, rng, 0); 
        
    uImp = piecewise(-intv<pt<intv, 1, 0); % critical to sit it in an interval...not accurate though
        fun_conT(uImp, rng, 22);    
end


%------------------------------------------------------------------------------------------ #1
if sel == 1
    
    %mag = 1/2;
    %ang = sym(pi);   %  -0.50  + j 0.000
    
    %mag = 1/2;
    %ang = -sym(pi);   %  -0.50  + j 0.00
    
    %mag = 1;
    %ang = sym(pi)/2;  % 0.00  + j 1.00
    
    %mag = 1;
    %ang = -sym(pi)/2;  % 0.00  - j 1.00
    
    %mag = 1;
    %ang = 5*sym(pi)/2;  % 0.00  + j 1.00
    
    %mag = sqrt(2);
    %ang = sym(pi)/4;  % 1.00  + j 1.00
    
    %mag = sqrt(2);
    %ang = 9*sym(pi)/4; % 1.00  + j 1.00
    
    cart = mag*( cos(ang) + 1j*sin(ang) );
    
    if imag(cart) < 0
        fprintf('--> cart:  %.2f  - j %.2f\n', real(cart), abs(imag(cart)));
    else
        fprintf('--> cart:  %.2f  + j %.2f\n', real(cart), imag(cart));
    end
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    %a = 5;
    %b = 0; % 5.00 exp( -j 0.00 )
    
    a = -2;
    b = 0; % 2.00 exp( j 3.14 )
    
    a = 0;
    b = -3; %  3.00 exp( -j 1.57 )
    
    a = 1/2;
    b= sqrt(3)/2; % 1.00 exp( j 1.05 )
    
    a = 1;
    b = 1; % 1.41 exp( j 0.79 )
    
    %z = a + 1j*b;
    
    z = (1-1j)^2; %2.00 exp( -j 1.57 )
    z = 1j*(1-1j);  % 1.41 exp( j 0.79 )
    z = (1+1j)/(1-1j); % 1.00 exp( j 1.57 )
    
    rads = angle(z);
    mag = abs(z);
    
    if rads > 0
        fprintf(' %.2f exp( j %.2f )\n', mag, rads);
    else
        fprintf(' %.2f exp( -j %.2f )\n', mag, abs(rads));
    end
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    %limit(exp(-2*pt), pt,inf) % 0
    %limit(heaviside(pt), pt, inf) % 1
    %limit( exp(1j*(2*pt+(pi/4))), pt, inf)
    syms T;
    syms N;
    
    % x(t) = exp(-2t) u(t)   just need 0--> inf  0 if t<0
    xt = exp(-2*pt);
    xt2 = xt^2;
    intg = int(xt2, pt);
    lmt = (1/(2*T))*( subs(intg, pt, T) - subs(intg, pt, 0) );
    pInf = limit(lmt, T, inf); 
    lmt = ( subs(intg, pt, T) - subs(intg, pt, 0) );
    eInf = limit(lmt, T, inf);
    %fprintf(' pInf:  %s     Einf:  %s  \n', pInf, eInf);
    
    xt = exp( 1j * (2*pt+(pi/4))  );
    xt2 = abs(xt); % key distinction
    intg = int(xt2, pt);
    lmt = (1/(2*T))*( subs(intg, pt, T) - subs(intg, pt, -T) );
    pInf = limit(lmt, T, inf);
    lmt =  subs(intg, pt, T) - subs(intg, pt, -T);
    Einf = limit(lmt, T, inf);
    %fprintf(' pInf:  %s     Einf:  %s  \n', simplify(pInf), simplify(Einf));
    
    xt = cos(pt);
    xt2 = xt^2;
    intg = int(xt2, pt);
    lmt = (1/(2*T))* ( subs(intg, pt, T) - subs(intg, pt, -T) );
    Pinf = limit(lmt, T, inf);
    lmt = ( subs(intg, pt, T) - subs(intg, pt, -T) );
    Einf = limit(lmt, T, inf);
    %fprintf('Pinf:  %s     Einf:  %s\n', simplify(Pinf), simplify(Einf));
    
    xt = (1/2)^i_n;
    xt2 = xt^2;
    temp = subs(xt2, i_n, N);
    sig = (1/(2*N+1))*symsum(temp, N, 0, N);
    Pinf = limit(sig, N, inf);
    Einf = symsum(xt2, i_n, 0, inf);
    %fprintf('Pinf:  %s     Einf:  %s\n', simplify(Pinf), simplify(Einf));
    
    xt = exp( 1j* (  (sym(pi)/(2*i_n)) + (sym(pi)/8) ) );
    xt2 = abs(xt)^2;
    temp = subs(xt2, i_n, N);
    sig = (1/(2*N+1))*symsum(temp, N);
    pInf = limit( sig, N, inf);
    Einf = symsum(xt2, i_n, -inf, inf);
    %fprintf('Pinf:  %s     Einf:  %s\n', simplify(Pinf), simplify(Einf)); wrong P = 1, E = inf  
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    rng = [-10,10];
    og = piecewise(i_n<-2, 0, i_n>4, 0, 1);
    
    trans = subs(og, i_n, -i_n);
    
    fun_disN2(og, rng, trans, rng, 0);
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    rng = [-10, 10];
    og = piecewise(pt<3, 0, 1);
    
    %trans = subs(og, pt, 1-pt);
    trans = subs(og, pt, 1-pt) + subs(og, pt, 2-pt);
    %trans = subs(og, pt, 3*pt);
    %trans = subs(og, pt, pt/3);
    fun_conT2(og, rng, trans, rng, 0);
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    rng = [-10,10];
    og = heaviside(i_n);
    comp = og - subs(og, i_n, i_n-4);
        %fun_disN(comp, rng,0);
        
    xt = sin(pt/2);
        %fun_conT(xt, rng, 0);
    
    trans = piecewise(i_n>=3,1,0)*(1/2)^i_n;
        fun_disN(trans, rng,0);
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    rng = [-3,3];
    og = piecewise(pt<-2, 0,...
                   -2<=pt<-1, pt+1,...
                   -1<=pt<0, 1,...
                   0<=pt<1, 2,...
                   1<=pt<2, 2-pt,...
                   pt>2, 0);
    ptlst = [ -2,0;
              -2,-1;
              -1,0;
              -1,1;
              0,1;
              0,2;
              1,2;
              1,1;
              2,0];
   %oga = subs(og, pt, pt-1);
   %ogb = subs(og, pt, -pt+2);
   %ogc = subs(og, pt, 2*pt+1);
   %ogd = subs(og, pt, 4-(pt/2));
   %oge = (og + subs(og, pt, -pt));
   %oge = (og + subs(og, pt, -pt))*heaviside(pt);
    oge1 = og * heaviside(pt);
    oge2 = subs(og, pt, -pt) * heaviside(pt);
   fun_conT2(oge1, rng, oge2, rng, 0);
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    rng = [-9,9];
    og = piecewise( in == -4, -1,...
                    in == -3, -1/2,...
                    in == -2, 1/2,...
                    -1 <= in <= 2, 1,...
                    in == 3, 1/2,...
                    0);
    %fun_disN(og, rng, 0);  %  og  is good
    oga = subs(og, in, in-4);
    ogb1 = subs(og, in, in+3);
    ogb2 = subs(og, in,3-in);
    ogc = subs(og, in, 3*in);
    ogd = subs(og, in, 3*in+1);
    oge = piecewise(-5<=in<=3, 1, 0) * og;
    fun_disN2(og, rng, oge, rng, 0);
end


%------------------------------------------------------------------------------------------ #49
if sel == 49
    syms x y z z1 z2
    x = sym(3);
    y = sym(4);
    %z = x + 1j*y;
    %z = ( sym(1) + sym(1j))^5;
    
    z1 = sym(sqrt(3))+sym(1j)^3;
    z2 = sym(1)-sym(1j);
    z = z1*z2;
    expand(z)
    r = abs(z)
    th = angle(z)
end


%------------------------------------------------------------------------------------------ #55
if sel == 55
    a = symsum( exp(1j*in*sym(pi)/2), in, 0, 9);
    b = symsum( exp(1j*in*sym(pi)/2), in, -2, 7);
    c = symsum( ((1/2)^in)*exp(1j*in*sym(pi)/2), in, 0, inf);
    d = symsum( ((1/2)^in)*exp(1j*in*sym(pi)/2), in, 2, inf)
    
end
    
    
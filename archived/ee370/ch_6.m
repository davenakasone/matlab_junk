%{
    chapter 6 : stability, routh-hurwitz...might want to use that TI

    1  :  trying out this guys function...ex6.1a
    2  :  trying out this guys function...ex6.1b
    3  :  still trying the function...pe6.1
    4  :  looking at special case, 0 in frist col, ex6.2
    5  :  ex6.8 using getting a closed loop tf   feedback(), pole()
    6  :  be able to get eigen()   ...they are the poles, pe6.4   use det(sI-A)==0
    7  :  hw, p3 basic Routh
    8  :  hw, p4 Roth with auxillary / all 0 in a row
    9  :  hw, p7 unity feedback ... tf() -> feedback() -> pole()
    10 :  hw, p24 finiding range of k
    11 :  hw, p26 find oscillation in fucked up diagram
    12 :  hw, p30

%}
close all;
clc;
select = 12;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    syms epsilon;
    my_den = (s + 2) * (s + 3) * (s + 5);
    pretty(expand(my_den))
    my_den_coeffs = coeffs(my_den);
    ra = f_routh(my_den_coeffs, epsilon) % seems ok
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    syms epsilon;
    my_den = s^3 + 10*s^2 + 31*s + 1030;
    pretty(expand(my_den))
    my_den_coeffs = coeffs(my_den);
    ra = f_routh(my_den_coeffs, epsilon) % seems ok
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    syms epsilon;
    my_den = 3*s^7 + 9*s^6 + 6*s^5 + 4*s^4 + 7*s^3 + 8*s^2 + 2*s + 6;
    pretty(expand(my_den))
    my_den_coeffs = coeffs(my_den);
    ra = f_routh(my_den_coeffs, epsilon) % seems bad, can't take lcd
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    syms epsilon;
    my_den = s^5 + 2*s^4 + 3*s^3 + 6*s^2 + 5*s +3;
    pretty(expand(my_den))
    my_den_coeffs = coeffs(my_den);
    rootz = roots(my_den_coeffs)
    ra = f_routh(my_den_coeffs, epsilon) % 
end


%------------------------------------------------------------------------------------------
if (select == 5)
    numg=128;
    deng=[1 3 10 24 48 96 128 192 0];
    G=tf(numg,deng);
    T=feedback(G,1)
    poles=pole(T)
end


%------------------------------------------------------------------------------------------
if (select == 6)
    A=[2, 1, 1; 1, 7, 1; -3, 4, -5]; % just need the A matrix...
Eig=eig(A)
end


%------------------------------------------------------------------------------------------
if (select == 7)
    H_num = [1, -2];
    H_den = [1, -2, 4, -3, 2, -3];
    Ts = tf(H_num, H_den);
    my_poles = pole(Ts)
    safety_check = roots(H_den)
end


%------------------------------------------------------------------------------------------
if (select == 8)
    H_num = [1, 2, 7, 21];
    H_den = [1, -2, 3, -6, 2, -4];
    Ts = tf(H_num, H_den);
    my_poles = pole(Ts)
    safety_check = roots(H_den)
end


%------------------------------------------------------------------------------------------
if (select == 9)
    syms s;
    getter = (s+2)*(s+3)*(s+4)*(s+5)
    pretty(expand(getter))
    pretty(expand(fliplr(coeffs(getter))))
    G_den = [1, 14, 71, 154, 120];
    G_num = 584;
    safety_check = (G_num/getter) / (1 + (G_num/getter));
    pretty(simplify(expand(safety_check)))
    G = tf(G_num, G_den);
    T = feedback(G,1)
    poles = pole(T)
    check_p = roots([1, 14, 71, 154, 704])
end


%------------------------------------------------------------------------------------------
if (select == 10)
    syms s;
    syms k;
    Gs = k * (s + 1) / ( s * (s + 2) * (s + 3) );
    Hs = (s + 5) / (s + 7);
    Ts = Gs / (1 + Gs * Hs);
    pretty(simplify(expand(Ts)));
end


%------------------------------------------------------------------------------------------
if (select == 11)
    syms s;
    syms K;
    G = s * (s + 3) * (s + 7);
    Ts = (K * G) / (1 + s*G + K*G);
    pretty(simplify(expand(Ts)))
end


%------------------------------------------------------------------------------------------
if (select == 12)
    syms s;
    syms k;
    Gs = k / (s * (s+1) * (s+2) * (s+6));
    Ts = Gs / (1 + Gs);
    pretty(simplify(expand(Ts)));
    Tsk = subs(Ts, k, 224/9);
    pretty(simplify(expand(Tsk)));
    my_tf = tf(224, [9, 81, 180, 108, 224]);
    my_poles = pole(my_tf)
    
end
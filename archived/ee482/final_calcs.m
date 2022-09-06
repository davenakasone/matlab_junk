%{
    ee482 final
    
    1  : 
    
   2
  
%}
close all;
clear;
clc;
select = 3;  %        <-----CHANGE


%------------------------------------------------------------------------------------------
if (select == 1)
    syms z;
    syms n;
    Hz = ((1 - 1/(2*z)) * (1 - 1/z)) / (1 + 1/(4*z) - 1/(8*z^2));
    [numz, denz] = numden(Hz);
    Hz_zeros = roots(sym2poly(numz))
    Hz_poles = roots(sym2poly(denz))
    [r,p,k] = residuez(sym2poly(numz), sym2poly(denz))
    zplane(sym2poly(numz), sym2poly(denz))
    ts = 1;
    tf_Hz = tf(sym2poly(numz), sym2poly(denz), ts);

    Hn = iztrans(Hz, n)
    check = simplify(ztrans(Hn, z) - Hz)
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms z;
    syms n;
    Hz = -1*((1/(3*z)) - 4/(3*z^2)) / (1 - 1/(4*z) + 1/(8*z^2));
    [numz, denz] = numden(Hz);
    Hz_zeros = roots(sym2poly(numz))
    Hz_poles = roots(sym2poly(denz))
    [r,p,k] = residuez(sym2poly(numz), sym2poly(denz))
    zplane(sym2poly(numz), sym2poly(denz))
    ts = 1;
    tf_Hz = tf(sym2poly(numz), sym2poly(denz), ts);

    Hn = iztrans(Hz, n)
    check = simplify(ztrans(Hn, z) - Hz)
end

%------------------------------------------------------------------------------------------
if (select == 3)
    x_1 = [2, 3, 2, 1];
    x_2 = [1, 1, 1, 1];
    h = conv(x_1, x_2)
    
end

%------------------------------------------------------------------------------------------
if (select == 99)
    
end
%%%%%%%%~~~~~~~~END>

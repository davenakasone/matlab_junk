%{
    chapter 9 : root locus, design

    1  :  ex9.1
    2  :  class ex 

    rlocus()
    rlocusplot()
    sgrid()
    step()
    Control System Designer

%}
format compact;
clear;
close all;
clc;
select = 2;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    syms K;
    zeta = 0.174;
    theta = pi - acos(zeta); % angle from positive x-axis
    fprintf("\nif zeta = %0.3f, then theata must be:  %0.3f rad  ,  %0.3f deg\n",...
        zeta, theta, rad2deg(theta));
    
    Gs_num = 1;
    Gsn = 1;
    Gs_den =  (s + 1) * (s + 2) * (s + 10);
    Gsd = sym2poly(Gs_den);
    Gs = Gs_num / Gs_den;
    KTs = K*Gs / (1 + K*Gs);
    %pretty(expand(Gs));
    tf_raw = tf(Gsn, Gsd);
    %f_rlocus(tf_raw);
    
    % stay on the zeta line...K = 164.6, -0.694 +/- j3.926
    KK = 164.6;
    kts = subs(KTs, K, KK); % try the new K
    [nn, dd] = numden(kts); 
    dd_vec = sym2poly(dd);
    rootz = roots(dd_vec); % update roots
    
    % looks good, now see the uncompensated error
    e_inf = 1/(1 + subs(KK*Gs, s, 0)); % using e(oo) = lim (1 / (1 + KGs)) s->0 ... 1 /Kp
    fprintf("\nerror uncompensated:  %0.3f\n", e_inf);
    
    
    
    
    comp_Hs_num = (s + 0.1);
    comp_Hs_den = s * Gs_den;
    comp_Hs = comp_Hs_num / comp_Hs_den;
    %pretty(expand(comp_Hs));
    compHsn = [1, 0.1];
    compHsd = [1, 13, 32, 20, 0];
    tf_comp = tf(compHsn, compHsd);
    %f_rlocus(tf_comp);
    
end


%------------------------------------------------------------------------------------------
if (select == 2)
    raw_tf = tf(1, [1, 2, 0]);
    f_rlocus(raw_tf);
end

%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
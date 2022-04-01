%{
    ee360 ch9

    #1      ex9.5       laplace trans
    #2      ex9.8       ilaplace trans
    #3      ex9.13      roc with linearity
    #4      p9.7        poles and possibilities
    #5      p9.9
    #6      p9.21a
    #7      p9.21b
    #8      p9.21c
    #9      p9.21d    whole transform looks off
    #10     p9.21g
    #11     p9.22b
    #12     p9.22c
    #13     p9.22d
    #14     p9.22f    what the fuck is the -t or t?
    #15     p9.22g    wtf -t ?
    #16     p9.26      ???
    #17     p9.27           ?
    #18     p9.28   just some pole-zero graphs
    #19     chw.1   determine transfer function
    #20     chw.2   find Hdb and graph it with meshc
    #21     chw.3   option 1 str8 Bode
    #22     chw.3   option 2 freq
    #23     chw.4   roots->poles/zeros  ->poly -> sys fun
    #24     chw.5   tf2zp -> zp2tf -> original
    #25     chw.6   zplane()
    #26     chw.7   residue
    #27     chw.8   residue to step response

    #28     9.29
    #29     9.30
    #31     9.31
    #33     9.33   good residue one


%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');


                sel = 31;  % CHANGE CHANGE CHANGE
                
                
global sig; syms sig; assume(sig, 'real'); % resistivity , sub in and out with const
global alpha; syms alpha; assume(alpha, 'real');
global beta; syms beta; assume(beta, 'real');
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var
global in; syms in; assume(in, {'real', 'integer'}); % index n
global s; syms s; % transform holder " j w0 "
global t; syms t; assume(t,'real'); % transform holder
global freq; syms freq; assume(freq, 'real'); % in Hz  f = 1 / T   or N
global omg; syms omg; assume(omg, 'real'); % omega , angular freq ... 2 pi f
global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global ik; syms ik; assume(ik, {'real', 'integer'}); % index k, usually for convolution
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, usually sub with inf
global T; syms T; assume(T, 'real');  % period bound, continous case
global t0; syms t0; assume(t0, 'real'); % arbitrary t0, usually offset
global tau; syms tau; assume(tau, 'real'); % dummy tau, usually for convolution


aaa = 1.5;% for chw
global a; syms a;
w00 = 5;% for chw
global w0; syms w0;

ee = cls_EE330_helper();
const = cls_CONST();
%publisher('ee330_hw4.m');
%------------------------------------------------------------------------------------------ #1
if sel == 1
    x1 = dirac(t);
    x2 = (-4/3) * exp(-t) * heaviside(t);
    x3 = (1/3) * exp(2*t) * heaviside(t);
    Xs = laplace(x1+x2+x3, t, s);
    pretty(Xs);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    Xs = 1/((s+1)*(s+2));
    xt = ilaplace(Xs,s,t);
    pretty(xt);  % not complete...can't tell unless you know original params...3 possibilties
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    Xs1 = 1/(s+1); % sig > -1
    Xs2 = 1/((s+1)*(s+2)); % sig > -1
    Xs = simplify(Xs1-Xs2)   %...rational cancel,  ROC now sig>-2  it grew
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    Xs = (s-1)/((s+2)*(s+3)*(s^2 + s + 1));
    num = 1;
    den = [-2,-3];
    
    temp_p = [1,1,1];  % handle s^2 + s + 1
    [rtz]= roots(temp_p);
    rtz = transpose(rtz);
    den = cat(2, den, rtz);
    
    %fun_zer_pol(num,den);
    
    pretty(ilaplace(Xs,s,t));
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    Xs = (2*(s+2))/(s^2+7*s+12);
    num=-2;
    den=[-3,-4];
    %fun_zer_pol(num,den);
    pretty(ilaplace(Xs,s,t))
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    xt = exp(-2*t)*heaviside(t)+exp(-3*t)*heaviside(t);
    Xs = laplace(xt,t,s);
    pretty(Xs);
    fun_zer_pol(-5/2,[-2,-3],1,-2)
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    xt = exp(-4*t)*heaviside(t)+exp(-5*t)*heaviside(t)*sin(5*t);
    Xs = laplace(xt,t,s);
    simplifyFraction(Xs)
    [Nu,De] = numden(Xs);
    zs = roots([1,15,70]);
    zs = transpose(zs);
    tem = roots([1,10,50]);
    ds = cat(2,-4,transpose(tem));
    fun_zer_pol(zs,ds,1,-4)
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    xt = exp(2*t)*heaviside(-t) + exp(3*t)*heaviside(-t);
    Xs = (-1)*laplace(xt,t,s);
    pretty(simplifyFraction(Xs));
    [Nu,De] = numden(Xs)
    z_cof = coeffs(Nu, s)
    p_cof = coeffs(De, s)
    z_rot = transpose(roots(z_cof))
    p_rot = transpose(roots(p_cof))
    %fun_zer_pol(5/2,[2,3],2,2);
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    
    Hs = ((s-2)^2  -  (s+2)^2)/((s+2)^2 * (s-2)^2);
    pretty(simplify(Hs));
    
    
    xt1 = t*exp(-2*t);
    xt1_int = int(xt1*exp(-s*t), t);
    xt1_ul = limit(xt1_int, t, inf);
    xt1_ll = limit(xt1_int, t, 0);
    Xs1 = xt1_ul - xt1_ll;
    
    xt2 = t*exp(2*t);
    xt2_int = int(xt2*exp(-s*t), t);
    xt2_ul = limit(xt2_int, t, 0);
    xt2_ll = limit(xt2_int, t, -inf);
    Xs2 = xt2_ul - xt2_ll;
    
    Xs = (1/(s+2)^2)-(1/(s-2)^2);
    %pretty(expand((s-2)^2 - (s+2)^2 ));
    %pretty(expand(((s+2)^2)*((s-2)^2)));
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    %fun_zer_pol(z_rot,p_rot,3,[-2,2]);
end


%------------------------------------------------------------------------------------------ #10
if sel == 10
    xt = heaviside(t)-heaviside(t-1);
    Xs = laplace(xt,t,s);
    fun_zer_pol(0,0,4,999);
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    Xs = s/(9+s^2);
    xt = ilaplace(Xs,s,t);
    
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,2,0);
end


%------------------------------------------------------------------------------------------ #12
if sel == 12
    Xs = (s+1)/(9+(s+1)^2);
    xt = ilaplace(Xs);
    
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,2,-1);
end


%------------------------------------------------------------------------------------------ #13
if sel == 13
    Xs = (s+2)/(s^2+7*s+12);
    pretty(partfrac(Xs));
    xt = ilaplace(Xs,s,t);
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,3,[-4,-3]);
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    Xs = (s+1)^2 / (s^2 - s + 1);
    %pretty(partfrac(Xs));
    xt = ilaplace(Xs,s,t);
    %pretty(xt);
    temp = ilaplace(xt,s,t);
    pretty(temp);
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,1,1/2);
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    Xs = (s^2 - s + 1)/(s+1)^2;
    %pretty(partfrac(Xs));
    xt = ilaplace(Xs,s,t);
    %pretty(xt);
    temp = ilaplace(xt,s,t);
    %pretty(temp);
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,1,-1);
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    Xs = ((1/(s+2))-(1/(s-3)))*exp(s);
    %pretty(partfrac(Xs));
    xt = ilaplace(Xs,s,t);
    %pretty(xt);
    temp = ilaplace(xt,s,t);
    %pretty(temp);
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,3,[-2,3]);
end


%------------------------------------------------------------------------------------------ #17
if sel == 17
    Xs = 16/(2+2*s+s^2);
    %pretty(partfrac(Xs));
    xt = ilaplace(Xs,s,t);
    %pretty(xt);
    temp = ilaplace(xt,s,t);
    %pretty(temp);
    [num,den]=numden(Xs);
    z_cof = coeffs(num,s,'All');
    p_cof = coeffs(den,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,1,-1);
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    num = 2;
    den = [-2,-1,1];
    %fun_zer_pol(num,den,2,-2);
    %fun_zer_pol(num,den,1,1);
    %fun_zer_pol(num,den,3,[-2,-1]);
    %fun_zer_pol(num,den,3,[-1,1]);
    fun_zer_pol(num,den,4,0);
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    pretty(Hs);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    pretty(Hss);
end


%------------------------------------------------------------------------------------------ #20
if sel == 20
    s = sig + 1j*omg;
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hs_mag = abs(Hs);
    %pretty(Hs_mag);
    Hss_mag = subs(Hs_mag, [a,w0], [aaa,w00]);
    %pretty(Hss_mag);
    
    bufS = 1;
    bufO = 1;
    bufH = 1;
    pts = 20;
    fil = 128;
    pos = [20,20,700,700];
    rng_sig = [-10,10];
    rng_omg = [-10,10];
    rng_Hmag = [0,0];
    tis = '20 log10 ( Hmag ) vs s = \sigma + j \omega';
    xax = linspace(rng_sig(1)-bufS, rng_sig(2)+bufS, fil);
    yax = linspace(rng_omg(1)-bufO, rng_omg(2)+bufO, fil);
    sig_ax = linspace(rng_sig(1)-bufS, rng_sig(2)+bufS, pts);
    omg_ax = linspace(rng_omg(1)-bufO, rng_omg(2)+bufO, pts);
    
    Hdb = zeros(pts,pts);
    [Sig,Omg] = meshgrid(sig_ax, omg_ax);
    for row = 1:pts
        for col = 1:pts
            Hdb(row,col)= 20*log10(subs(Hss_mag, [sig,omg], [Sig(row,col),Omg(row,col)]));
            if Hdb(row,col)< rng_Hmag(1)
                rng_Hmag(1) = Hdb(row,col);
            end
            if Hdb(row,col) > rng_Hmag(2)
                rng_Hmag(2) = Hdb(row,col);
            end
        end
    end
    zax = linspace(rng_Hmag(1)-bufH, rng_Hmag(2)+bufH, fil);
    
    figure('Name', 'question 2',...
           'Position', pos,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    view(125,30); % CHANGE
    title(tis, 'FontSize', 16);
    xlabel('real(s) = \sigma' , 'FontSize', 14);
    ylabel('imag(s) = \omega' , 'FontSize', 14); 
    zlabel('Hdb = 20log10(|H|', 'FontSize', 14);
    xlim([rng_sig(1)-bufS, rng_sig(2)+bufS]);
    ylim([rng_omg(1)-bufO, rng_omg(2)+bufO]);
    zlim([rng_Hmag(1)-bufH, rng_Hmag(2)+bufH]);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rng_sig(2)+bufS , 0          , 0                ,'y.', 'markersize', 20, 'linewidth', 10);  
    plot3(0           , rng_omg(2)+bufO, 0                ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0           , 0              , rng_Hmag(2)+bufH ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_sig(2)+bufS, 0              , 0               , '+ \sigma', 'FontSize', 14);
    text(0              , rng_omg(2)+bufO, 0               , '+ \omega', 'FontSize', 14);
    text(0              , 0              , rng_Hmag(2)+bufH, '+ |H| dB', 'FontSize', 14);
    
    plot3(-3/2, 0, 0,'bo', 'LineWidth', 3);
    text(-2, 0, 0, 'zero (-1.5,0)', 'FontSize', 10);
    plot3(-3/2, 5, 0,'rx', 'LineWidth', 3);
    text(-4, 6, 0, 'pole (-1.5,5)', 'FontSize', 10);
    plot3(-3/2, -5, 0,'rx', 'LineWidth', 3);
    text(-4, -6, 0, 'pole (-1.5,-5)', 'FontSize', 10);
    for ct = 1:fil
        plot3(-3/2,0,zax(ct),'b.', 'MarkerSize', 1);
        plot3(-3/2,5,zax(ct),'b.', 'MarkerSize', 1);
        plot3(-3/2,-5,zax(ct),'b.', 'MarkerSize', 1);
    end
    
    colormap jet;
    colorbar;
    hand = meshc(Sig,Omg,Hdb, 'FaceAlpha', .3, 'LineWidth', 1, 'FaceColor', 'texturemap');
    hand(2).LineWidth = 2;
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    [num,den] = numden(Hss);
    
    z_cof = double(coeffs(num,s,'All'));
    p_cof = double(coeffs(den,s,'All'));
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    fun_zer_pol(z_rot,p_rot,1,-3/2);
    
    H = tf(z_cof, p_cof);
    %bode(H);
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    pos = [20,20,700,700];
    pts = 200;
    rng_w = [-10,10];
    buf_w = 1;
    rng_mag = [0,0];
    buf_m = 1;
    rng_phs = [0,0];
    buf_p = 1;
    
    w = linspace(rng_w(1),rng_w(2),pts);
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
    H = freqs(b_cof, a_cof, w);
    
    Hmag = zeros(1,pts);
    Hphs = zeros(1,pts);
    for ct = 1:pts
        Hmag(1,ct) = 20*log10(abs(H(ct)));
        Hphs(1,ct) = (180/pi)*angle(H(ct));
        if Hmag(1,ct) < rng_mag(1) 
            rng_mag(1) = Hmag(1,ct);
        end
        if Hmag(1,ct) > rng_mag(2)
            rng_mag(2) = Hmag(1,ct);
        end
        if Hphs(1,ct) < rng_phs(1) 
            rng_phs(1) = Hphs(1,ct);
        end
        if Hphs(1,ct) > rng_phs(2)
            rng_phs(2) = Hphs(1,ct);
        end
    end
    ax_w = linspace(rng_w(1)-buf_w, rng_w(2)+buf_w, pts);
    ax_mag = linspace(rng_mag(1)-buf_m, rng_mag(2)+buf_m, pts);
    ax_phs = linspace(rng_phs(1)-buf_p, rng_phs(2)+buf_p, pts);
    figure('Position', pos);
    
    subplot(2,1,1);
    hold on;
    grid on;
    title('MAGNITUDE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('H dB = 20log10(|H|');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_mag(1)-buf_m, rng_mag(2)+buf_m]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_mag, ax_mag  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_mag(2)+buf_m, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hmag,'r-', 'LineWidth', 2);
    
    subplot(2,1,2);
    hold on;
    grid on;
    title('PHASE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('phase in °');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_phs(1)-buf_p, rng_phs(2)+buf_p]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_phs, ax_phs  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_phs(2)+buf_p, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hphs,'r-', 'LineWidth', 2);
end


%------------------------------------------------------------------------------------------ #23
if sel == 23
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    b_rot = roots(b_cof)
    b_pol = poly(b_rot);
    a_cof = double(coeffs(den,s,'All'));
    a_rot = roots(a_cof)
    a_pol = poly(a_rot);
    hand = zplane(b_pol, a_pol, 'g');
    hand.LineWidth = 2;
    %Hs = b_poly/a_poly
    %{
    z_cof = double(coeffs(num,s,'All'));
    p_cof = double(coeffs(den,s,'All'));
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    H = tf(z_cof, p_cof)
    %fun_zer_pol(z_rot,p_rot,1,-3/2);
    %}
end


%------------------------------------------------------------------------------------------ #24
if sel == 24
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
        H = tf(b_cof,a_cof)
        [zer,pol,con] = tf2zp(b_cof,a_cof)
        [b_out, a_out] =zp2tf(zer,pol,con)
        H = tf(b_out,a_out)
end


%------------------------------------------------------------------------------------------ #25
if sel == 25
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
    zplane(b_cof,a_cof);
end


%------------------------------------------------------------------------------------------ #26
if sel == 26
    Hs = (s+a)/((a^2 + w0^2)+2*s*a+s^2);
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    pretty(expand(Hss));
    H_lap = ilaplace(Hss,s,t);
    fprintf('the transfer function in the time domain:\n');
    pretty(H_lap);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
    
    pts = 100;
    pos = [20,20,700,700];
    rng_lap = [0,0];
    rng_res = [0,0];
    rng_tim = [0,5];
    buf_lap = 1;
    buf_res = 1;
    buf_tim = 1;
    
    ax_tim = linspace(rng_tim(1)-buf_tim, rng_tim(2)+buf_tim, pts);
    time = linspace(rng_tim(1), rng_tim(2), pts);
    [R, P, K] = residue(b_cof, a_cof);
    Hres = R.'*exp(kron(P,time));
    rng_res = [ min(Hres,[],'All'), max(Hres,[],'All')];
    ax_res = linspace(rng_res(1)-buf_res, rng_res(2)+buf_res, pts);
    Hlap = double(subs(H_lap, t, time));
    rng_lap = [ min(Hlap,[],'All'), max(Hlap,[],'All')];
    ax_lap = linspace(rng_lap(1)-buf_lap,rng_lap(2)+buf_lap,pts);
    
    figure('Position', pos);
    hold on;
    grid on;
    title('residue (blue)   Laplace(red)', 'FontSize', 26);
    view(2); 
    xlabel('time (t)', 'FontSize', 16);
    ylabel('h(t)', 'FontSize', 16);   
    xlim([rng_tim(1)-buf_tim, rng_tim(2)+buf_tim]);
    ylim([rng_res(1)-buf_res, rng_res(2)+buf_res]);
    plot(ax_tim      , 0 * ax_tim, 'k', 'linewidth', 1);
    plot(0 * ax_res, ax_res  , 'k', 'linewidth', 1);        
    plot(rng_tim(2)+buf_tim, 0      , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0      , rng_res(2)+buf_res, 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(time,Hlap,'r--', 'LineWidth', 4);
    plot(time,Hres,'b-', 'LineWidth', 1.5);
end

%------------------------------------------------------------------------------------------ #27
if sel == 27
    Hs = (s+a)/(s*((a^2 + w0^2)+2*s*a+s^2));
    Hss = subs(Hs, [a,w0], [aaa,w00]);
    
    [num,den]=numden(Hss);
    z_cof = double(coeffs(num,s,'All'));
    p_cof = double(coeffs(den,s,'All'));
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    H = tf(z_cof, p_cof);
    fun_zer_pol(z_rot,p_rot,1,0);
    
    pretty(expand(Hss));
    H_lap = ilaplace(Hss,s,t);
    fprintf('the transfer function in the time domain:\n');
    pretty(H_lap);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
    
    
    pts = 100;
    pos = [20,20,700,700];
    rng_lap = [0,0];
    rng_res = [0,0];
    rng_tim = [0,5];
    buf_lap = 1;
    buf_res = 1;
    buf_tim = 1;
    
    ax_tim = linspace(rng_tim(1)-buf_tim, rng_tim(2)+buf_tim, pts);
    time = linspace(rng_tim(1), rng_tim(2), pts);
    [R, P, K] = residue(b_cof, a_cof);
    Hres = R.'*exp(kron(P,time));
    rng_res = [ min(Hres,[],'All'), max(Hres,[],'All')];
    ax_res = linspace(rng_res(1)-buf_res, rng_res(2)+buf_res, pts);
    Hlap = double(subs(H_lap, t, time));
    rng_lap = [ min(Hlap,[],'All'), max(Hlap,[],'All')];
    ax_lap = linspace(rng_lap(1)-buf_lap,rng_lap(2)+buf_lap,pts);
    
    figure('Position', pos);
    hold on;
    grid on;
    title('residue (blue)   Laplace(red)', 'FontSize', 26);
    view(2); 
    xlabel('time (t)', 'FontSize', 16);
    ylabel('h(t)', 'FontSize', 16);   
    xlim([rng_tim(1)-buf_tim, rng_tim(2)+buf_tim]);
    ylim([rng_res(1)-buf_res, rng_res(2)+buf_res]);
    plot(ax_tim      , 0 * ax_tim, 'k', 'linewidth', 1);
    plot(0 * ax_res, ax_res  , 'k', 'linewidth', 1);        
    plot(rng_tim(2)+buf_tim, 0      , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0      , rng_res(2)+buf_res, 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(time,Hlap,'r--', 'LineWidth', 4);
    plot(time,Hres,'b-', 'LineWidth', 1.5);
end


%------------------------------------------------------------------------------------------ #28
if sel == 28
    xt = exp(-t)*heaviside(t);
    ht = exp(-2*t)*heaviside(t);
    xs = laplace(xt,t,s); %   1/(s + 1)    roc > -1
    hs = laplace(ht,t,s); % 1/(s + 2)   roc > -2
    
    ys = xs*hs;
    [num, den] = numden(ys);
end


%------------------------------------------------------------------------------------------ #29
if sel == 29
    x1t = heaviside(t);
    x1s = laplace(x1t,t,s);
    y1t = (1-exp(-t)-t*exp(-t))*heaviside(t);
    y1s = laplace(y1t,t,s);
    hs = simplify((1/x1s)*y1s);
    
    yt = (2-3*exp(-t)*exp(-3*t))*heaviside(t);
    ys = laplace(yt,t,s);
    %simplify(ys)
    xs = ys/hs;
    xt = ilaplace(xs,s,t);
    
    tst = (6*s+6)/(3*s+s^2);
    tstt = ilaplace(tst,s,t)
end


%------------------------------------------------------------------------------------------ #31
if sel == 31
    hs = 1/(s^2 - s -2);
    [hnum,hden] = numden(hs);
    
    z_cof = coeffs(hnum,s,'All');
    p_cof = coeffs(hden,s,'All');
    z_rot = transpose(roots(z_cof));
    p_rot = transpose(roots(p_cof));
    
    fun_zer_pol(-1,1/2, 999, 999);
    %fun_zer_pol(double(z_rot),double(p_rot), 3, [-1,2]);
    %fun_zer_pol(double(z_rot),double(p_rot), 1, 2);
    %fun_zer_pol(double(z_rot),double(p_rot), 2, -1);
end


%------------------------------------------------------------------------------------------ #33
if sel == 33
    hs = -2/((s^2 + 2*s +2)*(s-1));
    [hnum, hden] = numden(hs);
    zcof = double(coeffs(hnum, s, 'All'));
    pcof = double(coeffs(hden, s, 'All'));
    simplifyFraction(hs)
end
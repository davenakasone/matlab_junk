%{
    ch9, ch10, close out ee330

    #1      9.1  ex    conducting bar sliding, still, both  ...use general
    #2      9.1  pp    force, current, other emf resultants
    #3      9.2a pp    loop rotates
    #4      9.2b pp    loop and field rotate
    #5      9.3  ex    use reluctance to find Vemf of coils on a core
    #6      9.4  pp    capacitor and displacement current
    #7      9.5  ex    just some complex number review
    #8      9.6  ex    turn a vector into phasor   just some    exp(1j*omg*t)
    #9      9.6  pp    more phasor and instantaneous form on vectors
    #10     9.7  ex    satisfy maxwell with phasors
    #11     10.1 pp    mechanics of a wave
    #12     10.2 pp    wave is in a general medium, decouple
    #13     hw, 10.3   basic wave molesting
    #14     hw, 10.7   basic waves
    #15     hw, 10.15  basic waves
    #16     hw, 10.22  losless medium
    #17     hw, 10.24  good dielectric
    #18     hw, 10.31  it is a good dielectric
    #19     hw, #2       ...
    #20     10.4 pp    infer info from wave
    #21     10.8 ex    time avg power
    #22     10.8 pp    time avg power, free space
    #23     hw, 10.57   wave and power
    #24     hw, 10.63   2 interfaces
    #25     ee320, hw13, 10.96
    #26     ee320, hw13, 10.100

%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
format shortE; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]


                        select = 26;  % CHANGE CHANGE CHANGE
                        
                        
% Universal constants                
global null; null = double.empty();  % nulled array of size 1x1
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global j; j = sqrt(-1); % just overload it here so no problems occur later

% Greeks
global alp; syms alp; assume(alp, 'real'); % real phase, real z
global bet; syms bet; assume(bet, 'real'); % imag phase, imag z
global chie; syms chie; assume(chie, 'real'); % polarization
global chim; syms chim; assume(chim, 'real'); % magnetization
global ep; syms ep; assume(ep, 'real'); % as in ep = ep0 * epr    or D = ep * E
global ep0; syms ep0; assume(ep0, 'real'); % sub with const.ep0   permittivity of free space
global eta; syms eta; % intrinisic impedance
global gam; syms gam; % phase coeff
global lam; syms lam; assume(lam, 'real'); % wave length, eigen val, 
global mu; syms mu; assume(mu, 'real'); % as in mu = mu0 * mur    or B = mu * H , friction coeff
global mu0; syms mu0; assume(mu0, 'real'); %  sub with const.mu0   permiability of free space
global sig; syms sig; assume(sig, 'real'); % resistivity , mean, 
global omg; syms omg; assume(omg, 'real'); % omega , angular freq ... 2 pi f
global phi; syms phi; assume(phi, 'real'); % phase constant
global tau; syms tau; assume(tau, 'real'); % convolution dummy, torque
global tht; syms tht; assume(tht, 'real'); % any angle theta

% Single Vars (specific)
global amp; syms amp; assume(amp, 'real'); % amplitude of any function
global freq; syms freq; assume(freq, 'real'); % in Hz  f = 1 / T  = 2 pi / omg
global phs; syms phs; assume(phs, 'real'); % phasor term, (omg*t - beta*rz), pi/2, ect
global c1; syms c1; % integration constant 1
global c2; syms c2; % integration constant 2
global c3; syms c3; % integration constant 3
global c4; syms c4; % integration constant 4
global pL; syms pL; assume(pL, 'real'); assume(pL >= 0); % rho_line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); assume(pS >= 0); % rho_surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); assume(pV >= 0); % rho_volume charge denisty in C / m^3

% Single Vars (general)
global k; syms k; assume(k, {'real', 'integer'}); % index k
global n; syms n; assume(n, {'real', 'integer'}); % index n, iztrans
global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, DT period
global r; syms r; assume(r, 'real'); %  z = r exp(j tht)
global s; syms s; % j omg, laplace
global t; syms t; assume(t,'real'); % time , ilaplace
global t0; syms t0; assume(t0, 'real'); % arbitrary t0 offset
global T; syms T; assume(T, 'real');  % CT period, integration limit
global U; syms U; assume(U, 'real');   % U( X, Y)  ...f(Zxy) = U(X,Y) + j V(X,Y)  --> real(f)     
global V; syms V; assume(V, 'real');   % V( X, Y)  ...f(Zxy) = U(X,Y) + j V(X,Y)  --> imag(f)
global X; syms X; assume(X, 'real'); % real Z
global Y; syms Y; assume(Y, 'real'); % imag Z
global Z; syms Z;   % do not compound, z-trans
global Zxy; syms Zxy; Zxy = X + 1j*Y; % compound at will

% Coords       phi --> fi
global rx; global ry; global rz; % "rectangular x"   "rectangular y"  "rectangular z"       
global cr; global cf; global cz; % "cylindrical rho" "cylindrical fi" "cylindrical z"       
global sr; global st; global sf; % "spherical radius" "spherical theta" "spherical fi"   
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); assume(cr >= 0);
syms cf; assume(cf, 'real'); assume(cf >= 0);
syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); assume(sr >= 0);
syms st; assume(st, 'real'); assume(st >= 0);
syms sf; assume(sf, 'real'); assume(sf >= 0);

global arx; syms arx; assume( arx, 'real'); % x_hat (rec)
global ary; syms ary; assume( ary, 'real'); % y_hat (rec)
global arz; syms arz; assume( arz, 'real'); % z_hat (rec)
global acr; syms acr; assume( acr, 'real'); % r_hat (cyn)
global acf; syms acf; assume( acf, 'real'); % phi_hat (cyn)
global acz; syms acz; assume( acz, 'real'); % z_hat (cyn)
global asr; syms asr; assume( asr, 'real'); % r_hat (sph)
global ast; syms ast; assume( ast, 'real'); % theta_hat (sph)
global asf; syms asf; assume( asf, 'real'); % phi_hat (sph)

global Arx; syms Arx; assume( Arx, 'real'); % x_hat comp (rec)
global Ary; syms Ary; assume( Ary, 'real'); % y_hat comp (rec)
global Arz; syms Arz; assume( Arz, 'real'); % z_hat comp (rec)
global Acr; syms Acr; assume( Acr, 'real'); assume(Acr >= 0); % r_hat comp (cyn)
global Acf; syms Acf; assume( Acf, 'real'); assume(Acf >= 0); % phi_hat comp (cyn)
global Acz; syms Acz; assume( Acz, 'real'); % z_hat comp (cyn)
global Asr; syms Asr; assume( Asr, 'real'); assume(Asr >= 0); % r_hat comp (sph)
global Ast; syms Ast; assume( Ast, 'real'); assume(Ast >= 0);% theta_hat comp (sph)
global Asf; syms Asf; assume( Asf, 'real'); assume(Asf >= 0);% phi_hat comp (sph)

ee = cls_EE330_helper();
const = cls_CONST();
mat = cls_math('Themistocles');
%ep0 = const.ep0;
%------------------------------------------------------------------------------------------ #1
if select == 1 
    % bar can slide freely over 2 rails...const B in the area
    
    %{
    %a) changing field, static loop
        bar = .08;
        B = [0,0,.004*cos(t*10^6)];
        intg = dot( diff(B,t,1), [0,0,1]);
        Vemf = -int(int(intg, rx, 0, .06), ry, 0, .08);
        fprintf('V_emf = %s  V\n', Vemf);
    %}
    
    %{
    %b)  motional emf...bar moves, field const
        u = [0,20,0];
        B = [0,0,.04];
        intg = dot(cross(u,B),[1,0,0]);
        Vemf = int(intg, rx, .06, 0); % take the current that opposes
        fprintf('V_emf = %.3f  V\n', Vemf);
    %}
    
    %
    %c) both...   take it in parts, or take the total flux and use Faraday's Law
        u = [0,20,0];
        B = [0,0,.004*cos(t*10^6 - ry)];
        % dy/dt = u...
        
        intg = dot(B,[0,0,1]);
        f1 = int(intg, rx, 0, .06);
        f2 = int(f1, ry, 0, ry);
        flx = (-1)*diff(f2, t, 1);
        % dy/dt = u(t)      dy = u(t) dt      y = t u
        flux = subs(flx,ry, t*u(2))  
    %}
end


%------------------------------------------------------------------------------------------ #2
if select == 2
    B = [0,0,.5];
    R = 20;
    l = .1;
    u = [8,0,0];
    
    intg = dot(B,[0,0,-1]);
    f = int(int(intg, ry, 0, l), rx, 0, rx);
    % dx/dt = u(t)    dx = u(t) dt   x = t|u| + C    x=ut
    fl = subs(f, rx, t*u(1));
    flu = diff(fl, t, 1);
    emf = -1*flu; %  .4 V
    I = emf/R; % .02 mA
    Fm = [-I*l*norm(B),0,0]; %   - 1 mN x\hat
    pow = I^2 *R;  %  .8 mw
end


%------------------------------------------------------------------------------------------ #3
if select == 3
    R = .1;
    ll_z = 0;
    ul_z = .03;
    ll_r = 0;
    ul_r = .04;
    B_0 = .05;
    B_rec = [0,B_0,0]; % transform the mag field
    B_cyn = [B_0*sin(cf), B_0*cos(cf),0];
    w = 100*pi;
    
    u = [0,w * cr,0]; % w cr acf
    u_cp_B = cross(u, B_cyn);
    dl = [0,0,1]; % dz acz
    u_cp_B_dp_dl = dot(u_cp_B, dl);
    intg = subs(u_cp_B_dp_dl, cr, ul_r); % -.2 sin(cf) dz
    flux = int(intg, cz, ll_z, ul_z);
    
    % if  d{cf}/dt = w(t)  cf = wt + C    , C=pi/2
    fluxx = subs(flux, cf, w*t+sym(pi)/2);
    t1 = .001; % 1ms
    Vemf1 = double(subs(fluxx, t, t1));
    fprintf('t= 1ms , Vemf = %.3f  mV\n', Vemf1*(1e3));
    i = fluxx/R;
    t2 = 3e-3;
    i3 = subs(i, t, t2);
    fprintf('t= 3 ms, i = %.3f A\n', i3);
end


%------------------------------------------------------------------------------------------ #4
if select == 4
    t1 = 1e-3; % 1ms
    t2 = 3e-3; % 3ms
    R = .1;
    ll_z = 0;
    ul_z = .03;
    ll_r = 0;
    ul_r = .04;
    B_0 = .02*t;
    B_rec = [B_0,0,0]; % transform the mag field to rec
    B_cyn = [B_0*cos(cf), -B_0*sin(cf),0];
    w = 100*pi;
    
    % take the motional part
    mot_dl = [0,0,1]; % dz acz
    mot_u = [0,ul_r*w,0]; % CD*w == CD * d{cf}/dt
    mot_u_cp_B = cross(mot_u, B_cyn);
    mot_intg = dot(mot_u_cp_B, mot_dl);
    mot_flx = int(mot_intg, cz, ul_z, ll_z);
    mot_flux = subs(mot_flx, cf, w*t+sym(pi)/2);
    
    % take the transformer part
    trn_dS = [0, 1, 0]; % dr dz acf
    B_dt = diff(B_cyn, t, 1);
    trn_B_dp_dS = dot(B_dt, trn_dS);
    trn_intg = trn_B_dp_dS;
    trn_int = int(int(trn_intg, cr, ll_r, ul_r), cz, ll_z, ul_z);
    trn_flx = (-1)*trn_int;
    trn_flux = subs(trn_flx, cf, w*t+sym(pi)/2);
    
    % t = 1ms
    Vemf1 = subs(mot_flux, t, t1) + subs(trn_flux, t, t1);
    fprintf('@ t = 1 ms, Vemf = %.3f  uV\n', Vemf1*(10e6));
    
    % t = 3ms
    Vemf2 = subs(mot_flux, t, t2) + subs(trn_flux, t, t2);
    i=Vemf2/R;
    fprintf('@ t = 3 ms, i = %.3f  mA\n', i*(1e3));
end


%------------------------------------------------------------------------------------------ #5
if select == 5
    sa = 4*(10^-4); % 4 cm^2
    N1 = 500;
    N2 = 300;
    V1 = 120;
    % there is not much you have to do... V1 = -N1*flx  , V2 = -N2*flx... V2 = N2/N1 * flx
    V2 = (N2/N1)*V1; % 72 V
end


%------------------------------------------------------------------------------------------ #6
if select == 6
    E = [0, 20*cos(omg*t-50*rx), 0];
    D = ep0 .* E;  
    J_d = diff(D, t, 1);
    H = [0,0, (-1)*int(J_d(2), rx) ]; % match cross product
    cp_E = ee.getCurlRec(E);
    H_dt = (-1)*mu0*diff(H,t,1);
    w = sqrt((5/2)*(1000/(const.ep0*const.mu0)));
    fprintf('omg = %.3f  Grad/sec\n', w/(10^9));
end


%------------------------------------------------------------------------------------------ #7
if select == 7
    a = (1j*conj(3-(1j)*4))/( (-1+1j*6) * (2+1j)^2 ); %  .16216 - j .027027
    a_mag = abs(a); % .1644
    a_ang = rad2deg(angle(a)); % -9.46°
    b  = sqrt( (1+1j)/(4-1j*8) ); % .2325 + j .3226
    b_mag = abs(b); % .3976
    b_ang = rad2deg(angle(b)); % 54.422°
    
    d = (1j)^3 * ( (1+1j)/(2-1j) )^2;
    d_re = real(d); % .24
    d_im = imag(d); % .32
    
    c = 6*exp(1j*pi/6) + 1j*5 -3 + exp(1j*pi/4);
    c_re = real(c); % ??? fuck this book
    c_im = imag(c); % 8.707
end


%------------------------------------------------------------------------------------------ #8
if select == 8
    w = 10^8;
    phs = sym(pi)/3;
    syms phz; assume(phz, 'real');
    
    Ap = [0,0, 10*cos(omg*t-10*rx+phz)];
    Bs = [20/(1j), 10*exp(1j*2*phz*rx), 0];
    A_re = real( [0,0,10*exp(1j*omg*t)*exp(1j*(phz-10*rx))] );
    %simplify(Ap-A_re); % [0,0,0] is good
    
    A_s = [0,0,10*exp(1j*(phz-10*rx))];
    %simplify(Ap-real(exp(1j*omg*t) .* A_s)); % [0,0,0] is good
    
    B = real( exp(1j*omg*t) .* Bs);
    simplify(B-real( exp(1j*omg*t) .* Bs)); % [0,0,0] is good
end


%------------------------------------------------------------------------------------------ #9
if select == 9
    % switch the sin of the P\vec
    phs1 = rx-3*sym(pi)/4; 
    w1 = 10;
    P = [0,2*cos(omg*t+phs),0];
    P_sub = subs(P, [omg, phs], [w1, phs1]);
    Ps = [0,2*exp(1j*phs1),0];
    Ps_sub = subs(Ps, [omg,phs], [w1,phs1]);
    simplify(P_sub - real( exp(1j*w1*t) * Ps_sub) ); % [0,0,0] means it worked
    
    Q_mag = sin(sym(pi)*ry);
    Q_phs = rx;
    Qs = [exp(1j*Q_phs)*Q_mag, exp(1j*Q_phs)*Q_mag, 0];
    Q = simplify(real(exp(1j*omg*t) .* Qs));
    simplify(Q-real( exp(1j*omg*t) .* Qs)); % [0,0,0] means good
end


%------------------------------------------------------------------------------------------ #10
if select == 10
    syms H0;
    w = 10^6;
    ph = beta * cz;
    E_mag = 50/cr;
    H_mag = H0/cr;
    expp = exp(1j*w*t);
    
    E = [ 0 , E_mag * cos(w*t+ph) , 0 ];
    Es = [ 0 , E_mag * exp(1j*ph) , 0 ];
    testE = simplify(E-real(expp.*Es));
    subs(testE, [t,beta], [666,999]) % good phasor vector
    
    H = [ H_mag * cos(w*t+ph),0,0];
    Hs = [ H_mag * exp(1j*ph), 0, 0];
    testH = simplify(H-real(expp.*Hs));
    subs(testH, [H0, t, beta, cz], [111,22,333,44]); % good phasor
end


%------------------------------------------------------------------------------------------ #11
if select == 11
    rxx = ones(1,100);
    amp = .1;
    omg = 2e8;
    u = const.c;
    bet = omg/u; %  2/3    rad/m
    lam = (2*pi)/bet;  %  9.425  m
    T = (2*pi)/omg; %   31.42 ns
    phs = omg * t - bet * rx;
    %t0 = T/8 - lam/(8*u); % =0
    t0 = T/8
    Hy = amp * cos(phs)
    Hyy = zeros(1,100);
    for ct = 1:100
        Hyy(ct)=amp*cos(omg*t0-bet*rxx(ct));
    end
    plot(rxx, Hyy)
    H = [ 0, Hy, 0];
end


%------------------------------------------------------------------------------------------ #12
if select == 12
    epr = 8;
    ep = const.ep0 * epr;
    mur = 2;
    mu = const.mu0 * mur;
    omg = 10^8;
    alp = 1/3;
    
    t1 = (alp/omg)^2 * (2/(ep*mu)) + 1;
    sig_omg_ep2 = t1^2 -1;
    bet = double(omg * sqrt( ((mu*ep)/2) * ( sqrt(1+sig_omg_ep2) + 1 ) )); % 1.374 rad/m
    gam = alp + j*bet;
    sig = double(   sqrt(((bet^2 + alp^2)/(omg*mu))^2   -   (omg^2 * ep^2))    );
    
    eta = double(  sqrt( (j*omg*mu) / (sig+j*omg*ep) )   );
    eta_th = angle(eta); %  .23794  rad
    eta_thd = rad2deg(eta_th); % or  13.63°
    eta_mag = abs(eta);  % 177.72
    th_loss = double(sqrt(sig_omg_ep2)); % .5154
    chk1 = atan(th_loss)/2; % .23794
    etaM = double(sqrt(mu/ep)/(1+(sig/(omg*ep))^2)^(1/4)); % just a check =  177.72
    etaA = double(atan((sig/(omg*ep)))/2); % just a check, =  .23794 rad
    
    u = omg/bet; % 7.276e7 m/s;
    H0 = 1/(2*eta_mag); % .002814
    Hy = (H0 * sin(omg*t-bet*rz-eta_th))
    H = [0, Hy, 0];
end


%------------------------------------------------------------------------------------------ #13
if select == 13
    H0 = .4;
    omg = 10^8;
    alp = 0;
    T = (2*pi)/omg;
    freq = 2*pi*omg;
    sig = 0;
    mu = const.mu0;
    ep = const.ep0;
    
    bet = omg*sqrt(mu*ep); % 1/3  rad / m
    lam = (2*pi)/bet; % 18.85 m
    u = omg/bet; % = c
    syms y;
    Hy = -H0 * cos( omt*t + bet*y);
    
    time = [ 2e-9, 3e-9, 4e-9, 10e-9];
end


%------------------------------------------------------------------------------------------ #14
if select == 14
    freq = 50 * 10^6;
    omg = 2 * pi * freq;
    mu = 2.1*const.mu0;
    ep = 3.6*const.ep0;
    sig = .08;
    E0 = 6;
    
    alp = double(const.alp_cal(omg, mu, ep, sig)); % 5.4106
    bet = double(const.bet_cal(omg, mu, ep, sig)); % 6.129
    gam = alp + j*bet; %  5.4106e+00 + 6.1290e+00i
    lam = (2*pi)/bet;  %  1.0251e+00  m
    u = omg/bet;  % 5.1257e+07  m/s
    
    eta = double(sqrt( (j*omg*mu)/(sig+j*omg*ep) )); %   7.6021e+01 + 6.7110e+01i
    eta_mag = abs(eta); % 1.0141e+02
    eta_ang = angle(eta); % 7.2322e-01
    eta_angd = rad2deg(eta_ang); %  4.1437e+01  °

    H0 = E0/eta_mag; %   5.9168e-02    
end


%------------------------------------------------------------------------------------------ #15
if select == 15
    f = 10^9;
    omg = 2*pi*f;
    sig = 1;
    ep = 4*const.ep0;
    mu = 9*const.mu0;
    
    alp = double(const.alp_cal(omg, mu, ep, sig));   %  1.6882e+02
    bet = double(const.bet_cal(omg, mu, ep, sig)); % 2.1046e+02
    gam = alp + j * bet; %  1.6882e+02 + 2.1046e+02i
    lam = (2*pi)/bet;  %  2.9855e-02
    u = omg/bet;
    
    eta = double(sqrt( (j*omg*mu)/(sig+j*omg*ep) )); %  2.0545e+02 + 1.6480e+02i
    eta_mag = abs(eta);  %  2.6338e+02
    eta_ang = angle(eta); %  6.7606e-01  rad 
    eta_angd = rad2deg(eta_ang); %  3.8736e+01  °
end


%------------------------------------------------------------------------------------------ #16
if select == 16
    f = 40 * 10^6;
    omg = 2*pi*f;  %  2.5133e+08
    mu = const.mu0;
    ep = 4.5*const.ep0;
    Ex = 8;
    Ey = -6;
    sig = 0;
    alp = 0;
    bet = double(omg*sqrt(ep*mu)); % 1.7772e+00
    lam = 2*pi/bet; %  3.5355e+00
    u = omg/bet; %  1.4142e+08
    
    eta = double(sqrt(mu/ep));  %  1.7772e+02
    Hy = Ex/eta;  % 4.5016e-02
    Hx = Ey/eta; % 3.3762e-02
end


%------------------------------------------------------------------------------------------ #17
if select == 17
    H0 = 20e-3;
    mu = const.mu0;
    ep = const.ep0;
    alp = 12;
    bet = 12;
    freq = 10^6;
    omg = 2*pi*freq;  % 6.2832e+06
    sig = double((bet^2)/(pi*freq*mu)); % 36.48
    
    eta = (alp + j*alp)/sig; %  3.2899e-01 + 3.2899e-01i
    eta_mag = abs(eta); %  4.6526e-01
    eta_rad = angle(eta);  %  7.8540e-01
    eta_deg = rad2deg(eta_rad);  %  45°
    
    E0 = H0/eta_mag;  %  9.3052e-03     4.2987e-02 ???  
end


%------------------------------------------------------------------------------------------ #18
if select == 18
    freq = 12*10^6;
    omg = 2*pi*freq;
    mu = const.mu0;
    eta_mag = 24.6;
    eta_rad = pi/4;
    eta = eta_mag * exp(j*eta_rad); % 1.7395e+01 + 1.7395e+01i
    sig = double((omg*mu)/eta_mag^2); %  1.5657e-01
    alp = double(sqrt(pi*freq*mu*sig)); %  2.7235e+00
    bet = alp;
    lam = double(2*pi/bet); % 2.3071e+00
    u = omg/bet
end


%------------------------------------------------------------------------------------------ #19
if select == 19
    E0 = 10;
    freq = 10^9;  %  1 GHz
    omg = 2*pi*freq; %    6.2832e+09    rad/sec
    bet = 100*pi;  %     3.1416e+02   Np/m
    lam = 2*pi/bet;  %  2.0000e-02  m
    alp = 1;
    u = double(omg/bet); %  2e7   m/s
    
    mu = const.mu0;  % it is not magnetic
    ep = double((bet^2 - alp^2)/(mu*omg^2)); %  1.9894e-09   F/m
    temp1 = ((bet^2 + alp^2)/(omg*mu))^2;
    temp2 = (omg*ep)^2;
    sig = double(sqrt(temp1-temp2)); %  7.9577e-02
    test_alp = double(const.alp_cal(omg, mu, ep, sig));
    test_bet = double(const.bet_cal(omg, mu, ep, sig));
    %fprintf('difference of alphas: %.4f\n', test_alp-alp);
    %fprintf('difference of betas : %.4f\n', test_bet-bet);
    
    eta = double(const.eta_cal(omg, mu, ep, sig));% 2.5132e+01 + 7.9999e-02i
    eta_mag = abs(eta); % 2.5133e+01
    eta_rad = angle(eta);
    eta_deg = rad2deg(eta_rad); %  1.8238e-01
    
    loss_tan = sig/(omg*ep); %  6.3663e-03
    chk_lt = tan(2*eta_rad); %   6.3663e-03
    chkk = omg*ep;   %   1.2500e+01
    
    skin = 1/alp; %  1
    
    H0=E0/eta_mag; %  3.9789e-01
    
    %{
    ep_eqn = -(ep^3)*(omg^4 * mu * 4 * bet^2) + (ep^2)*(omg^2 * 4 * bet^4 - omg^4 * mu^2) - ep*(omg^2 * mu * 4 * alp^2)-4*alp^4;
    epz = double(solve(ep_eqn==0, ep));
    ep1 = epz(1) % no neg
    ep2 = epz(2) % seems legit
    ep3 = epz(3) % no imag
    
    ep = ep2;
    sig = double(omg*ep*sqrt( ((2*alp^2)/(omg^2 * mu * ep))^2 -1 ))
    alp_chk = double(const.alp_cal(omg, mu, ep, sig))   % 
    bet_chk = double(const.bet_cal(omg, mu, ep, sig))   % 
    %}
end


%------------------------------------------------------------------------------------------ #20
if select == 20
    % wave travels in ary direction , the medium is lossy
    ep = 4 * const.ep0;
    mu = const.mu0;
    sig = 10^-2;
    omg =pi*10^9;
    E0 = 30;
    %E = [0,0, E0 * cos(omg*t+(pi/4))];  @ t=0
    alp = double(const.alp_cal(omg, mu, ep, sig));    %   9.4153e-01   Np/m
    bet = double(const.bet_cal(omg, mu, ep, sig));     %  2.0965e+01   rad/m
    gam = alp + j*bet;
    E = [0,0, E0 * exp(-alp*ry) * cos(omg*t-bet*ry+(pi/4)) ];
    E12 = double(subs(E, [t, ry], [ (2e-9), 1] ));  %   2.7887e+00     V/m
    lam = 2*pi/bet;   %   2.9970e-01  m
    u = omg/bet; %  1.4985e+08  m/s
    
    % dist traveled in 10°   travels lam in 360...
    d10 = lam*(1/36);   %   8.3249e-03
    
    %amplitude reduced by 40%  
    y40 = log(.6)/(-alp); %   5.4255e-01  m
    
    % decouple for magnetic field   ...must be attenuating on -arx direction
    eta = double(const.eta_cal(omg, mu, ep, sig));
    eta_mag = abs(eta); %  1.8812e+02
    eta_rad = angle(eta); %   4.4879e-02   rad
    eta_deg = rad2deg(eta_rad); %  2.5714e+00  °
    H0 = E0/eta_mag; %   1.5948e-01
    H = [ -H0*exp(-alp*ry)*cos(omg*t-bet*ry+(pi/4)-eta_rad),0,0]; % E leads, so subtract
    H22 = double(subs(H, [ry,t], [2,2e-9])); %   2.2798e-02   arx  A/m
end


%------------------------------------------------------------------------------------------ #21
if select == 21
    % the medium is non-magnetic
    omg = 2*pi*10^7;
    bet = .8;
    alp = 0;
    gam = alp+j*bet;
    phs = -pi/2;
    E0 = 4;
    Ez = E0 * cos(omg*t-bet*rx+phs); %  shifted off the sin(wt-bx)
    E = [0,0,Ez]; % prop in arx, z-polarized   V/m
    
    % you can see the medium is not in free space, but probably lossless
    mu = const.mu0; % because non-magnetic   ...there is still some magnetisim
    ep = double(((bet*const.c)/omg)^2 *const.ep0);
    eta = double(sqrt(mu/ep)); %  9.8696e+01   ohm
    
    % find the mag field...it must be -ary
    H0 = -E0/eta; % -4.0528e-02
    H = [ 0, -H0*cos(omg*t-bet*rx+phs), 0]; % same phase
    pointing = cross(E,H);
    % [-(45631004489637*cos((4*rx)/5 - (4216574282663131*t)/67108864 + pi/2)^2)/281474976710656, 0, 0]
    T = 2*pi/omg; %  1.0000e-07
    pavg = simplify((1/T)*int(pointing, t, 0, T));
    pavgg = E0^2 / (2*eta);  %  8.1057e-02
    pavggg = [double((1/T)*int(pavgg, t, 0, T)),0,0];  %   8.1057e-02  watts   arx
    
    % the plane turned into a gradient.....2x+y=5    ->  arn = [2,1,0] / sqrt(5)
    arn = (1/sqrt(5)) .* [2,1,0];
    intg = dot(pavggg, arn);
    Pavg = intg * 100*(1/10^4);  % 725 uW
end


%------------------------------------------------------------------------------------------ #22
if select == 22
    H0 = .2;
    alp = 0;
    mu = const.mu0;
    ep = const.ep0;
    eta = sqrt(mu/ep);
    E0 = H0 * eta;
    pvg = [E0^2 / (2*eta),0,0];
    arn = (1/sqrt(2)) .* [1,1,0];
    pavg = dot(pvg,arn);
    Pavg10 = double(pavg*100*(1/10^4));  %  5.3315e-02  W
    area  = pi*(.05)^2;
    arnn = (area) .* [1,0,0];
    Pavg5 = double(dot(pvg, arnn));
end


%------------------------------------------------------------------------------------------ #23
if select == 23
    H0 = 4;
    alp = 0;
    bet = 5;
    sig = 0;
    phs = -pi/2;
    mu = const.mu0;
    ep = 4*const.ep0;
    omg = double(bet/sqrt(mu*ep)); %  7.5e8  rad/m
    u = omg/bet; %  1.5e8  m/s
    lam = 2*pi/bet; %  1.2566  m
    eta = double(sqrt(mu/ep));  %   1.8850e+02
    E0 = H0*eta; %  7.5398e+02
    
    et1 = sqrt(const.mu0/const.ep0);
    et2 = sqrt(mu/ep);
    Gam = (et2-et1)/(et2+et1); %  -1/3;
    Er0 = double(Gam*et1*H0);
    
    tau = (2*et2)/(et2+et1); %  2/3
    Et = double(tau * H0*et1); %  1.0053e+03
    Ht = double(Et/et2)
    pony = [Ht*Et,0,0];  %  5.3617e+03
end


%------------------------------------------------------------------------------------------ #24
if select == 24
    ep1 = const.ep0;
    mu1 = const.mu0;
    alp1 = 0;
    sig1 = 0;
    bet1 = 3;
    gam1 = alp1+j*bet1;
    omg1 = bet1/sqrt(ep1*mu1); %  9e8  rad/m
    u1 = omg1/bet1;  %  c
    lam1 = 2*pi/bet1; %  2.094 m
    eta1 = double(sqrt(mu1/ep1)); % 377
    E01 = 10;
    H01 = E01/eta1;  % .02653
    
    ep2 = 80*const.ep0;
    mu2 = const.mu0;
    sig2 = 4;
    lt2 = sig2/(omg1*ep2); %   2 pi 
    alp2 = double(const.alp_cal(omg1, mu2, ep2, sig2)); %  4.3936e+01
    bet2 = double(const.bet_cal(omg1, mu2, ep2, sig2)); %  5.1482e+01
    eta2 = double(const.eta_cal(omg1, mu2, ep2, sig2)); % 1.2711e+01 + 1.0848e+01i
    eta2_mag = abs(eta2); % 1.6710e+01
    eta2_rad = angle(eta2); %  7.0648e-01
    eta2_deg = rad2deg(eta2_rad); %   4.0478e+01  °
    
    Gam = (eta2-eta1)/(eta1+eta2);
    Gam_mag = abs(Gam); %   9.3482e-01
    Gam_deg = rad2deg(angle(Gam)); %  1.7670e+02  °
    
    tau = (2*eta2)/(eta2+eta1); %  6.6730e-02 + 5.3814e-02i
    tau_mag = abs(tau); %   8.5725e-02
    tau_deg = rad2deg(angle(tau)); %   3.8884e+01  °
end


%------------------------------------------------------------------------------------------ #25
if select == 25
    g_m = 5e-3;
    R_sig = 100e3;
    R_G1 = 47e6;
    R_G2 = 10e6;
    R_G = const.para(R_G1,R_G2); %  8.2456e+06
    Cc1 = .01e-6;
    R_S = 2e3;
    C_S = 10e-6;
    R_D = 4.7e3;
    R_L = 10e3;
    Cc2 = 1e-6;
    
    A_M = (-R_G/(R_G+R_sig))*g_m*const.para(R_D,R_L); %  -1.5795e+01      "   -15.80      "
    
    wp1 = 1/(Cc1*(R_sig+R_G)); %  1.1982e+01
    fp1 = wp1/(2*pi); %  1.9070e+00                                       "   1.91 Hz     "
    
    wp2 = (g_m+(1/R_S))/C_S; % 5.5000e+02
    fp2 = wp2/(2*pi); %  8.7535e+01                                       "   87.54  Hz   "
    
    wp3 = 1/(Cc2*(R_D+R_L)); % 6.8027e+01
    fp3 = wp3/(2*pi); %  1.0827e+01                                       "   10.83  Hz   "
    
    wz = 1/(C_S*R_S); % 50
    fz = wz/(2*pi); % 7.9577e+00                                          "   7.96  Hz    "
    
    fl = sqrt( fp1^2 + fp2^2 + fp3^2 -2*fz^2); %  8.7502e+01              "   87.50  Hz   "
end


%------------------------------------------------------------------------------------------ #26
if select == 26
    R_sig = 5e3;
    R_B1 = 33e3;
    R_B2 = 22e3;
    R_B = const.para(R_B1, R_B2); % 13200 " 13.2 k"
    R_E = 3.9e3;
    R_C = 4.7e3;
    R_L = 5.6e3;
    V_CC = 5;
    I_E = .3e-3;
    beta = 120;
    alph = beta/(beta+1);
    C_C1 = 1e-6;
    C_C2 = 1e-6;
    C_E = 20e-6;
    %r_e = const.sil_VTroom/I_E; % 8.3333e+01
    %r_e=25; %???
    
    % neglect the Early effects, solve bias:
    g_m = (I_E*beta)/((beta+1)*const.sil_VTroom); %  1.1901e-02
    r_pi = beta/g_m; %  1.0083e+04
    t1 = r_pi + (1+beta)*R_E;
    R_in = const.para(R_B, t1); % 1.2848e+04                      "  12.848 k ohm "
    r_e = r_pi/(1+beta); % 8.3333e+01
    
    % use the pi model to find the mid-band gain:
    t2 = (R_in/(R_in+R_sig));
    t3 = (r_pi/(r_pi+(1+beta)*R_E));
    t4 = const.para(R_C,R_L);
    A_M = -g_m*t2*t3*t4; %   -4.5798e-01                           "   -.458  V/V  "
    
    % solve the resistances
    R_C1 = const.para(R_B, r_pi)+R_sig; % 1.0717e+04
    t5 = const.para(R_B, R_sig);
    t6 = r_e + (t5/(beta+1));
    R_CE = const.para(R_E, t6); %  1.1010e+02
    R_C2 = R_C + R_L; % 10.3 k
    
    % solve time consts
    tau_C1 = R_C1*C_C1; % 1.0717e-02        " 10.72 ms "
    tau_CE = R_CE*C_E;  % 2.2021e-03        " 2.20  ms "
    tau_C2 = R_C2*C_C2; % 1.0300e-02        " 10.30 ms "
    
    fl = (1/(2*pi))*( (1/tau_C1) + (1/tau_CE) + (1/tau_C2) ); % 1.0258e+02  "  102.58 Hz "
    
end
    
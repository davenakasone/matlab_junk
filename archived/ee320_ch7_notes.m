%{
    ee320 ch7 transistor amplifiers

    #1      pp7.1   voltage amplifier
    #2      ex7.1   voltage gains nmos
    #3      ex7.2   bjt voltage amp
    #4      pp7.3   extension of ex7.2
    #5      pp7.4   NMOS small signal
    #6      pp7.5   NMOS  with Early
    #7      pp7.6   NMOS, work backwards
    #8      pp7.7   PMOS, NMOS, align g_m
    #9      pp7.8   PMOS
    #10     pp7.10  small signal model NMOS
    #11     pp7.15  BJT as amp
    #12     pp7.16  BJT analysis as voltage amplifier
    #13     pp7.19  BJT(pnp) with caps
    #14     pp7.20  BJT(npn) bias current
    #15     pp7.21  CS amp
    #16     ex7.8   CE amp
    #17     pp7.25  CG amp
    #18     pp7.26  CB amp
    #19     pp7.27  CB amp fed by coaxial cable
    #20     pp7.28  CD "source follower"
    #21     pp7.29  CD "source follower" as final stage in cascade
    #22     pp7.30  CC "emitter follower"
    #23     pp7.31 "darlington"  CC emitter followers paired
    #24     ex7.11  solve these bias arrangements simotaneously
    #25     pp7.32  looking at ex7.11, but with the bad Vgs bias...much more change
    #26     pp7.33  bias with 2 power sources
    #27     pp7.34  feedback from R_G
    #28     pp7.36  BJT classical w/ 2 power supplies
    #29     pp7.37  BJT classical w/ feedback

    #30     q7.85 HW1
    #31     q9.95 HW2
    #32     q7.120 HW3
    #33     q7.130 HW4
    #34     q7.134 HW5

%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
format shortE; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]

const = cls_CONST();
global t; syms t; assume(t,'real'); % transform holder
global w; syms w; assume(w, 'real'); % lazy omega
V_T = const.sil_VTroom;

                        select = 34;  % CHANGE CHANGE CHANGE

                        
%------------------------------------------------------------------------------------------ #1
if select == 1
    vdd = 1.8;
    rd = 17.5e3; 
    vtn = .4;
    kn = 4e-3;
    lam = 0; % ignore early effect
    
    % pt "A" is where NMOS goes cut-off to saturation   Vds < Vtn  it is cut off
    % Vds >= Vtn  in saturation.... Vds = Vt  on the border
    vgsa = vtn; % .4 V
    id = 0;
    vdsa = vdd - id*rd; % 1.8 V
    fprintf('vtc vgs to vds @ A = ( %.3f V, %.3f V )\n', vgsa, vdsa);
    
    % get the "B" ends the saturation, NMOS is going into triode
    vgsb = vtn + ( sqrt(2*kn*rd*vdd+1) - 1 ) / (kn*rd);
    vovb = vgsb - vtn;
    vdsb = vovb;
    fprintf('vtc vgs to vds @ B = ( %.3f V, %.3f V )\n', vgsb, vdsb);
    
    % at point "C" you need vgsc and vdsc, but infer from id
    vgsc = vdd; % because it is independent
    syms vdsc;
    Vdsc = vdd - kn * ( (vgsc-vtn)*vdsc - (1/2)*vdsc^2 )*rd - vdsc;
    slv = double(solve(Vdsc==0,vdsc));
    vdsc = slv(1); % .183 mV
    fprintf('vdsc @ C = %.3f mV\n', vdsc*(1e3));
end


%------------------------------------------------------------------------------------------ #2
if select == 2
    Vtn = .4;
    knp = .4e-3;
    WL = 10;
    lam = 0;
    kn = WL*knp;
    Vdd = 1.8;
    Rd = 17.5e3;
    Vgs = .6; % dc bias
    
    % for vgs = 0 ,  vds = 0, find Vov, Id, Vds, Av
    Vov = Vgs - Vtn; % .2 V
    Id = (1/2)*kn*Vov^2; %  .08 mA
    Vds = Vdd - Id*Rd; % .4 V
    
    % max allowable swing is amplitude of sinusiod v_gs
    Av = -1*kn*Vov*Rd; % -14 V/V      Vds > Vov   indeed saturation
    % anything over .2 V in magnitude can cause cut-off
    v_ds_max = .2;
    v_gs_max = v_ds_max/abs(Av); % 14.3 mV
    % v_gs_max << Vov  so operation is about linear
    
    % to remain in saturation at negative peak,   v_DS.min >= v_GS.max - V_tn
    v_gs_max = .2/(abs(Av)+1); %  13.3 mV   slightly lower
end


%------------------------------------------------------------------------------------------ #3
if select == 3
    I_S = 10^-15;
    R_C = 6.8e3;
    V_CC = 10;
    V_CE_edge = .3;
    
    % find bias voltage V_BE required to operate @ V_CE = 3.2 V   and I_C there
    V_CE = 3.2;
    I_C = (V_CC - V_CE)/R_C; % 1 mA
    syms vbe;
    eqn_V_BE = I_S * exp(vbe/const.sil_VTroom);
    V_BE = double(solve(eqn_V_BE == I_C, vbe)); % 690.8 mV
    
    % find A_v (voltage gain) at this point if sinusoid of amplitude 5 mV is on V_BE
    % assum linear operation
    A_v = (-1)*(V_CC - V_CE)/const.sil_VTroom; %  -272 V/V
    amp = 5e-3;
    v_ce_max = abs(A_v) * amp; % 1.36 V
    
    % find positive increment in v_BE (above V_BE) the transistor is put on edge of saturation
    I_C_edge = (V_CC - V_CE_edge) / R_C; % 1.43 mA
    % to increase I_C from 1 mA to 1.43 mA, solve delta_v_BE = V_T ln(1.43/1)
    incremnt = const.sil_VTroom * log(I_C_edge/I_C); % about 9 mV
    
    % find negative increment in v_BE that drives transitor to within 1% of cut-off 
    %(v_CE = .99 V_CC)
    I_C_cut = (V_CC - .99*V_CC)/R_C; % 14.7 uA
    % to decrease I_C from 1mA to 14.7 uA, delta_v_BE = V_T ln(.0147/1)
    inc = const.sil_VTroom * log(I_C_cut/I_C); % -105.5 mV
end


%------------------------------------------------------------------------------------------ #4
if select == 4
    I_S = 10^-15;
    I_C = 1e-3;
    V_CC = 10;
    V_CE_edge = .3;
    
    % find R_C so A_v is -320 V/V   just use the eqn:  A_v = -Ic Rc / VT
    A_v = -320;
    R_C = A_v * const.sil_VTroom / -I_C; % need 8k ohms
    
    % min v_CE , should not go below .3 V
    v_CE_min = (A_v * const.sil_VTroom) + V_CC; % 2 V     V_CE is 2 V @ Q
    % in order to go from 2 V to .3 V  
    v_CE_neg = v_CE_min - V_CE_edge; % 1.7 V
    % solve from implied agin
    inp_neg = v_CE_neg/abs(A_v); % 5.3 mV 
end


%------------------------------------------------------------------------------------------ #5
if select == 5
    V_DD = 5;
    R_D = 10e3;
    V_tn = 1;
    k_np = 20e-6;
    WL = 20;
    V_GS = 2;
    lam = 0;
    
    % infer
    k_n = k_np * WL;
    V_OV = V_GS - V_tn;
    
    % find I_D
    I_D = (1/2) * k_n * V_OV^2; % .2 mA
    
    % find V_DS
    V_DS = V_DD - I_D * R_D; % 3 V
    
    % find g_m
    g_m = 2 * I_D / V_OV; % .4 mA / V  
    
    % find voltage gain A_v 
    A_v = -k_n * V_OV * R_D; % -4 V/V
    A_v = -g_m * R_D; %  -4 V/V
    
    % v_gs = .2 sin(wt)    use small signal, min and max v_DS
    v_gs_mag = .2;
    v_gs = v_gs_mag * sin(w*t);
    v_ds = A_v * v_gs; % -.8 sin(t*w)  V
    v_ds_mag = abs(subs(v_ds, [t,w],[1,sym(pi)/2])); % .8
    
    % combine ac, dc for min/max values
    v_DS_min = double(V_DS - v_ds_mag); % 2.2 V
    v_DS_max = double(V_DS + v_ds_mag); % 3.8 V
    
    % components of i_D by trig identities   sin(wt)^2  =  1/2 - 1/2 cos(2wt)
    i_D = (1/2) * k_n * (V_GS - V_tn)^2 + k_n * (V_GS - V_tn) * v_gs + (1/2)* k_n * (v_gs)^2;
    pretty(rewrite(i_D, 'cos')); % the second harmonic is that cos(2wt) term
    fun_harm = (1/12500); % 8 uA
    sec_harm = (1/2)*(1/125000); % 4 uA
    distor = sec_harm / fun_harm % 5 percent second harmoic distortion
end


%------------------------------------------------------------------------------------------ #6
if select == 6
    k_np = 60e-6;
    WL = 40;
    k_n = k_np * WL;
    V_tn = 1;
    V_A = 15;
    lam = 1/V_A;
    
    % V_GS = 1.5 V  , find g_m  r_o
    V_GS = 1.5;
    g_m = k_n * (V_GS - V_tn); %  1.2  mA/V
    I_D = (g_m)^2 / (2 * k_n);% .3 mA
    r_o = abs(V_A) / I_D; %  50 k ohm
    
    % I_D = .5 mA  , now find g_m  r_o
    I_D = .5e-3;
    g_m = sqrt(2 * k_n * I_D); %  1.5492  mA / V
    r_o = abs(V_A) / I_D; % 30 k ohm
end


%------------------------------------------------------------------------------------------ #7
if select == 7
    I_D = .1e-3;
    g_m = 1e-3;
    k_np = 500e-6;
    
    % infer W/L  from regular g_m eqn...
    WL = ( g_m^2 ) / ( 2 * k_np * I_D ); % 10
    k_n = WL * k_np;
    V_OV = sqrt ( (2 * I_D) / k_n ); % .2 V
end


%------------------------------------------------------------------------------------------ #8
if select == 8
    rat = 1/.4; % 2.5   ...everything else is the same
end


%------------------------------------------------------------------------------------------ #9
if select == 9
    V_tp = -1;
    k_pp = 60e-6;
    WL = 16/.8;
    k_p = k_pp * WL;
    V_SG = 1.6;
    
    % find I_D
    I_D = (1/2) * k_p * ( V_SG - abs(V_tp) )^2;  %  216 uA
    g_m = sqrt( 2 * k_p * I_D);  % .72 mA/V
    g_m = (2*I_D)/(V_SG - abs(V_tp)); % same
end


%------------------------------------------------------------------------------------------ #10
if select == 10
    R_D = 10e3;
    R_G = 10e6;
    V_DD = 5;
    V_tn = .7;
    k_n = 1e-3;
    R_in = .5e6;
    g_m = -25;
    lam = 0; % neglect Early here
    
    I_D = ( g_m^2 ) / ( 2 * k_n )
end


%------------------------------------------------------------------------------------------ #11
if select == 11
    beta = 100;
    I_C = 1e-3;
    g_m = I_C / const.sil_VTroom; % 40 mA/V
    r_e = 1 / g_m;  % 25 ohm
    r_pi = (1 + beta) * r_e; % about 25 k ohm
end


%------------------------------------------------------------------------------------------ #12
if select == 12
    I_C = 1e-3;
    V_CC = 15;
    R_C = 10e3;
    beta = 100;
    v_be_mag = 5e-3;
    v_be = v_be_mag * sin(w*t);
    V_T = const.sil_VTroom;
    
    % find gain A_v = v_ce / v_be
    g_m = I_C / V_T; % 400
    A_v = (-I_C * R_C)/ V_T; % -400 
    v_ce_mag = abs(A_v * v_be_mag); % 2
    v_ce = v_ce_mag * sin(w*t);
    v_C_t = V_CC - I_C * R_C - v_ce; % 5 - 2*sin(t*w)
    
    % base current
    i_B_t = (I_C/beta) + ((I_C/V_T)*v_be)/beta;
end


%------------------------------------------------------------------------------------------ #13
if select == 13
    R_C = 7.5e3; % only real change from ex7.7
    R_E = 10e3;
    V_EE = 10;
    V_CC = -10;
    V_E = .7; % v_EB = .7 , V_EB = V_E - V_B ...base is grounded 
    beta = 100;
    alpha = beta/(1+beta);
    lam = 0; % ignore Early
    max_mag = 10e-3; % small signal requires amplitudes below 10 mV
    
    I_E = (V_EE - V_E)/R_E; % .93 mA
    I_C = alpha * I_E; % .92 mA
    V_C = V_CC + I_C * R_C; % -3.1 V
    
    g_m = I_C / const.sil_VTroom; % .0368  A/V
    r_e = const.sil_VTroom / I_E; % 26.88 ohm
    r_pi = beta/g_m; % 2.72 k ohm
    
    % use the T-model...it has an easy grounded base
    syms v_i;
    i_e = -v_i/r_e;
    v_o = -alpha*i_e*R_C;
    A_v = double(v_o / v_i); % 276.24 V/V
    
    max_vo = max_mag * A_v; % 2.76 V
end


%------------------------------------------------------------------------------------------ #14
if select == 14
    R_C = 8e3;
    V_CC = 10;
    R_B = 10e3;
    I = 1e-3; 
    beta = 100;
    alpha = beta/(1+beta);
    V_A = 100;
    v_BE = .7;
    
    % neglect early, find DC @ V_B, V_C, V_E   not sure how this is possible
    I_E = I;
    I_C = alpha * I_E; % .99 mA
    V_C = V_CC - I_C*R_C; % 2.1 V
    I_B = I_C/beta; % 9.9 uA
    V_B = -I_B*R_B; % -.1 V
    V_E = V_B - v_BE; % -.8 V
    
    g_m = I_C / const.sil_VTroom; % .04 A/V
    r_pi = const.sil_VTroom/I_B; %  2.525 k ohm
    r_e = const.sil_VTroom/I_E; % 25 ohm
    r_o = abs(V_A)/I_C; % 101 k ohm
    
    % z is grounded, x on v_sig, y on rest...use hybrid-pi model..small signal equiv amp
    syms v_sig;
    syms v_pi;
    R_sig = 2e3;
    R_L = 8e3;
    % inf cap means that caps can be short at any freq
    % I should be open when making small signal
    para = const.para(const.para(r_o, R_L),R_C);
    v_pi = v_sig/2; % KCL @ base(X)
    v_y = -g_m*v_pi*para;
    A_v_ro = double(v_y/v_sig); % -76.2 V/V
    % now neglect early (r_o)
    para = const.para(R_L,R_C);
    v_y = -g_m*v_pi*para;
    A_v = double(v_y/v_sig); % -79.2 V/V
    
    error = 100*((A_v - A_v_ro)/A_v_ro); % 3.9%
end


%------------------------------------------------------------------------------------------ #15
if select == 15
    I_D = .25e-3; % for bias on the MOSFET
    V_OV = .25;
    R_D = 20e3;
    R_sig = 100e3;
    R_L = 20e3;
    
    % R_in = oo  fyi...
    g_m = 2*I_D/V_OV;
    A_vo = -g_m*R_D; % -40 V/V
    R_o = R_D; % 20 k ohm
    A_v = -g_m * const.para(R_D, R_L); % -20 V/V
    G_v = A_v; % -20 V/V   ...see table
    
    % peak is limited to 10% of 2V_OV
    v_sig_max = 2*V_OV*.1; %  50 mV
    v_o_max = abs(v_sig_max * A_v); % 1 V
end


%------------------------------------------------------------------------------------------ #16
if select == 16
    beta = 100;
    alpha = beta/(beta+1);
    I_C = 1e-3;
    R_C = 5e3;
    R_sig = 5e3;
    R_L = 5e3;
    
    g_m = I_C/V_T;
    r_pi = alpha * V_T/I_C;
    r_e = alpha/g_m;
    R_in = (beta+1)*r_e; % 2.5 k ohm
    R_o = R_C; % 5 k ohm
    A_vo = -g_m*R_C; % -200 V/V
    
    A_v = -g_m*const.para(R_C,R_L); % -100 V/V
    G_v = -beta*( const.para(R_C,R_L) / ( R_sig + R_in ) ); % -33.33 V/V
    
    v_pi_mag = 5e-3;
    v_sig_mag = ( (R_in + R_sig) / R_in ) * v_pi_mag; % 15 mV
    v_o_mag = abs(v_sig_mag * G_v); % .5 V
end


%------------------------------------------------------------------------------------------ #17
if select == 17
    R_sig = 100;
    V_OV = .2;
    R_D = 2e3;
    % match R_sig...
    R_in = R_sig;
    I_D = V_OV/(2*R_sig); % 1 mA
    g_m = (2*I_D)/V_OV;
    G_v = R_D/(R_sig+(1/g_m)); % 10 V/V 
end


%------------------------------------------------------------------------------------------ #18
if select == 18
    I_C = 1e-3;
    R_C = 5e3;
    R_sig = 5e3;
    R_L = 5e3;
    
    g_m = I_C / V_T;
    r_e = V_T / I_C; % alpha about 1 ?
    
    R_in = r_e; % 25 ohm       ...yes, just make alpha=1 if unsure
    A_vo = g_m*R_C; %  200 V/V
    R_o = R_C;   %   5 k ohm
    
    A_v = g_m * const.para(R_C, R_L); %  100 V/V    ...notice these are positive gains
    G_v = const.para(R_C,R_L)/(R_sig + r_e); % .5 V/V   ...pathetic gain
end


%------------------------------------------------------------------------------------------ #19
if select == 19
    R_sig = 50;
    r_e = R_sig;
    I_C = V_T/R_sig; %  .5 mA   matched
    G_v = 40;
    R_tot = G_v * (R_sig + r_e); % 4 k ohm
end


%------------------------------------------------------------------------------------------ #20
if select == 20
    V_OV = .25;
    R_o = 100;
    g_m = 1/R_o;
    R_L = 1e3;
    v_sig = 1;
    v_i = v_sig;
    R_sig = 1e6;
    I_D = g_m*V_OV/2; % 1.25 mA
    v_o = v_i * ( R_L/(R_L + R_o)); % .91 V
    v_gs = v_i*(R_o/(R_o+R_L)); % .91 mV   ...very slick voltage division
end


%------------------------------------------------------------------------------------------ #21
if select == 21
    R_o = 200;
    k_np = .4e-3;
    V_OV = .25;
    g_m = 1/R_o;
    WL = g_m/(k_np*V_OV); % 50
    k_n = WL * k_np;
    I_D = g_m*V_OV/2; % .625 mA
    R_Lmin = 1e3;
    R_Lmax = 10e3;
    G_v_min = R_Lmin/(R_Lmin+R_o); %  .833 V/V
    G_v_max = R_Lmax/(R_Lmax+R_o); % .98 V/V
end


%------------------------------------------------------------------------------------------ #22
if select == 22
    beta = 100;
    alpha = beta/(beta+1);
    I_C = 5e-3;
    R_sig = 10e3;
    R_L = 1e3;
    
    g_m = I_C/V_T;
    r_e = alpha*V_T/I_C;
    r_pi = beta*V_T/I_C;
    A_vo = 1; % given
    G_vo = 1; % given
    R_o = r_e;
    R_in = (1+beta)*(r_e+R_L); % 101.5 k ohm
    R_out = r_e + (R_sig/(beta+1)); % 104 ohm
    A_v = R_L/(R_L+r_e);
    G_v = R_L/(R_L+r_e+(R_sig/(beta+1))); % .906 V/V
    
    v_pi_mag = 5e-3; % peak amplitude across junction
    v_o = R_L*v_pi_mag/r_e; % 1.01
    v_sig = v_o * ( (R_L+R_out)/R_L); % 1.12
end


%------------------------------------------------------------------------------------------ #23
if select == 23
    beta = 100;
    beta1 = beta;
    beta2 = beta;
    alpha = beta/(beta+1);
    I_E2 = 5e-3;
    I_C = ((1+beta)/beta)*I_E2;
    
    R_E = 1e3;
    R_sig = 100e3;
    
    r_e2 = V_T/I_E2;
    r_e1 = r_e2;
    
    R_in = (beta1 + 1)*(r_e1+(beta2+1)*(r_e2+R_E)); % 10.3 M ohm
    t1 = const.para(R_E, r_e2);
    t2 = (r_e1+(R_sig/(beta1+1)) )/(beta2+1);
    R_out = t1+t2; % 15 ohm
    t3 = (r_e1+(R_sig/(beta1+1)))/(beta2+1);
    G_v = R_E/(R_E+r_e2+t3); % .99
end


%------------------------------------------------------------------------------------------ #24
if select == 24
    I_D1 = .5e-3;
    V_tn1 = 1;
    k_n = 1e-3;
    V_DD = 15;
    V_SS = 0;
    R_G1 = 8e6;
    R_G2 = 7e6;
    lam = 0;
    V_D = 10; % pick the 2/3 and 1/3 drops
    V_S1 = 5;
    
    R_D = (V_DD-V_D)/I_D1; % 10 k ohm
    R_S = (V_S1 - V_SS)/I_D1; % 10 k ohm
    V_OV = sqrt(2*I_D1/k_n); % 1 V
    V_GS1 = V_tn1 + V_OV; % 2 V
    V_G1 = V_S1 + V_GS1; % 7 V
    check = V_DD * (R_G2/(R_G1+R_G2)); % 7 V   good voltage divider feeds gate
    
    % asses the change
    V_tn2 = 1.5;
    syms I_D2;
    syms V_GS2;
    I_D2 = (1/2)*k_n*(V_GS2 - V_tn2)^2;
    % second eqn is  V_G = V_GS2 + I_D2 R_S   ....7= V_GS2 + 10 ( I_D2 )
    equa = V_GS2 + R_S * I_D2;
    sols = double(solve(equa==7, V_GS2));
    pos_IDa = double(subs(I_D2, V_GS2, sols(1)));% .665 mA   not valid
    pos_IDb = double(subs(I_D2, V_GS2, sols(2)));% .455 mA
    
    change = 100*(pos_IDb - I_D1)/I_D1; % decreases 9%
end


%------------------------------------------------------------------------------------------ #25
if select == 25
    I_D1 = .5e-3;
    V_tn1 = 1;
    k_n = 1e-3;
    lam=0;
    
    V_GS = V_tn1 + sqrt(2*I_D1/k_n); % 2 V
    V_tn2 = 1.5; % V_GS will be const
    I_D2 = (1/2)*k_n*(V_GS-V_tn2)^2;
    
    change = (I_D2-I_D1)/I_D1; % it dropped 75%
end


%------------------------------------------------------------------------------------------ #26
if select == 26
    I_D = .5e-3;
    V_D = 2;
    V_tn = 1;
    k_n = 1e-3;
    V_DD = 5;
    V_SS = -5;
    
    V_GS = V_tn + sqrt(2*I_D/k_n); % 2 V
    V_S = -V_GS; % because gate is grounded
    R_D = (V_DD-V_D)/I_D; % 6 k ohm
    R_S = (V_S - V_SS)/I_D; % 6 k ohm
end


%------------------------------------------------------------------------------------------ #27
if select == 27
    I_D = .5e-3;
    V_DD = 5;
    k_n = 1e-3;
    V_tn = 1;
    lam = 0;
    
    V_GS = sqrt(2*I_D/k_n)+V_tn; % 2 V
    R_D = (V_DD-V_GS)/I_D; % 6 k ohm
    R_D_5 = R_D * 1.05;
    I_D2 = I_D*(R_D/R_D_5);
    
    V_GS2 = sqrt(2*I_D2/k_n)+V_tn;
    I_D3 = (V_DD - V_GS2)/R_D_5
    % no more iterations needed
end


%------------------------------------------------------------------------------------------ #28
if select == 28
    beta = 100;
    alpha = beta/(beta+1)
    I_E = 1e-3;
    swing_min = -2;
    swing_max = 2;
    V_CC = 10;
    V_EE = -5;
    V_B = 0;
    V_BE = .7;
    V_E = V_B - V_BE;
    R_B = 0; % only needed for coupled signal
    
    R_E = (V_E-V_EE)/I_E; % 4.3 k ohm   KVL the junction
    I_C = alpha*I_E;
    V_C = V_E + 2;
    R_C = (V_CC - V_C)/I_C;
end


%------------------------------------------------------------------------------------------ #29
if select == 29
    I_E = 1e-3;
    swing = 2; % +/-
    V_CE = 2.3;
    V_CC = 10;
    beta = 100;
    alpha = beta/(beta+1);
    
    V_E = 0;
    V_C = V_CE - V_E;
    V_BE = .7;
    V_B = V_BE - V_E;
    I_B = I_E/(beta+1);
    R_B = (V_C-V_BE)/I_B; %  161 k ohm
    R_C = (V_CC - V_BE-I_B*R_B)/I_E; % 7.7 k ohm
end


%------------------------------------------------------------------------------------------ #30
if select == 30   % hw1,7_85
    I_C = 5e-3;
    beta = 200;
    alpha = beta/(beta+1);
    R_sig = 10e3;
    R_L = 200;
    G_vo = 1;
    
    g_m = I_C / V_T; %  .2
    r_e = (V_T/I_C)* alpha; % 4.975
    r_pi = (V_T/I_C)*beta; % 1 k ohm
    
    R_in = (1+beta)*(r_e+R_L);  %  41.2 k ohm
    vb_vsig = R_in/(R_in+R_sig); % .805 V/V
    G_v = R_L/(R_L+r_e+(R_sig/(beta+1))); % .785 V/V
    G_vv = vb_vsig * (R_L/(R_L+r_e)); % same
    
    vbe_max = 10e-3;
    v_o_max = (vbe_max/r_e)*R_L; % .402 V
    v_sig_max = v_o_max/G_v; % .512 V
    
    R_out = r_e + R_sig/(beta+1); % 54.73 ohms
    G_vvv = G_vo * (R_L/(R_L+R_out)); % matches
    R_LL=150;
    G_vvvv = G_vo * (R_LL/(R_LL+R_out)); % .733
end


%------------------------------------------------------------------------------------------ #31
if select == 31   % hw2,7_95
    V_G = 4;
    R_S = 15e3;
    V_tn = .7;
    k_n1 = 5e-3;
    k_n2 = 1.5*k_n1;
    
    syms I_D;
    V_GS = 4 - I_D * R_S;
    eqn = (1/2)*k_n1*( V_GS - V_tn)^2;
    ids = double(solve(eqn == I_D, I_D)); %  .2011 or .241 mA
    eqnn = (1/2)*k_n2*( V_GS - V_tn)^2;
    idss = double(solve(eqnn == I_D, I_D)); %  .204 or .237 mA
end


%------------------------------------------------------------------------------------------ #32
if select == 32   % hw3,7_120
    R_in2 = 50;
    g_m2 = 1/R_in2; % .02
    g_m1 = g_m2;
    v_i_mag = 5e-3;
    R_cab = 50;
    
    i_D = v_i_mag * g_m1; % .1 mA
    v_d1_mag = i_D * R_cab; % 5 mV
    v_o_mag = 1;
    R_D = v_o_mag / i_D; % 10 k ohm
end


%------------------------------------------------------------------------------------------ #33
if select == 33   % hw4,7_130
    V_CC = 15;
    beta = 100;
    alpha = beta/(beta+1);
    R_C = 6.8e3;
    R_1 = 100e3;
    R_2 = 47e3;
    R_E = 3.9e3;
    R_L = 2e3;
    R_sig = 5e3;
    V_BE = .7;
    
    %
    V_B = V_CC * (R_2/(R_1+R_2)); % 4.8 V
    R_B = const.para(R_1,R_2); % 31.973 k ohm
    I_E = (V_B - V_BE)/(R_E+(R_B/(1+beta))); % .971 mA
    I_C = alpha * I_E;
    V_C = V_CC - I_C * R_C; % 8.46 V
    
    g_m = I_C / V_T; % .0385
    r_pi = (V_T/I_C)*beta; % 2.6 k ohm
    r_e = (V_T/I_C)*alpha; % 25.74 ohm
    R_in1 = const.para(R_B,r_pi); % 2.4 k ohm
    R_in2 = R_in1;
    vb1_vsig = (R_in1/(R_in1+R_sig)); %  .325
    v_b2_vb1 = -g_m*const.para(R_in2, R_C); % -68.33 V/V
    vo_vb2 = -g_m*const.para(R_C,R_L); % -59.46 V/V
    G_v = vb1_vsig*v_b2_vb1*vo_vb2; %  1319 V/V
    %}
    
    %{ 
    syms vb ib;
    vb = R_E * (1+beta) * ib + V_BE;
    lhs = (R_2/R_1)*(V_CC-vb)-ib*R_2;
    I_B = double(solve(lhs==vb, ib)); % 9.62 uA
    I_E = (beta+1)*I_B; % .971 mA
    I_C = alpha * I_E; % .962 mA
    V_C = V_CC - I_C * R_C; % 8.46 V
    V_B = (R_1*R_2/(R_1+R_2))*((V_CC/R_1)-I_B); % 4.49 V
    %}
end


%------------------------------------------------------------------------------------------ #34
if select == 34   % hw5,7_134
    V_CC = 3;
    R_sig = 10e3;
    R_B = 100e3;
    R_E = 1e3;
    betaa = 80;
    alphaa = betaa/(betaa+1);
    betab = 180;
    alphab = betab/(betab+1);
    V_BE = .7;
    
    I_Ea = (V_CC - V_BE)/ (R_E + R_B/(betaa+1)); % 1.03 mA
    V_Ea = I_Ea*R_E; % 1.03 V
    V_Ba = V_Ea + V_BE; % 1.73 V
    I_Ca = alphaa*I_Ea; % 1.02 mA
    g_ma = I_Ca/V_T; % .04
    r_ea = V_T/I_Ea; % 24.3 ohm
    r_pia = betaa*(V_T/I_Ca); % 1967 ohm
    
    I_Eb = (V_CC - V_BE)/ (R_E + R_B/(betab+1)); % 1.48 mA
    V_Eb = I_Eb*R_E; % 1.48 V
    V_Bb = V_Eb + V_BE; % 2.18 V
    I_Cb = alphab*I_Eb; % 1.47 mA
    g_mb = I_Cb/V_T; % .06
    r_eb = alphab/g_ma; % 24.5 ohm
    r_pib = betab/g_mb; % 3.054 k
    
    R_e = const.para(R_E,R_E); % 500 ohm
    R=R_e;
    R_ina = const.para(R_B,(1+betaa)*(r_ea+R_e)); % 29.8 k ohm
    R_inb = const.para(R_B,(1+betab)*(r_eb+R_e)); % 48.7 k ohm
    G_va = (R_ina/(R_ina+R_sig))*(R/(R+r_ea));% .714 V/V
    G_vb = (R_inb/(R_inb+R_sig))*(R/(R+r_eb));% .791 V/V
end
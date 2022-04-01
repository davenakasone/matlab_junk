%{
    1  : problem 1
    2  : problem 2
    3  : problem 3
    4  : problem 4
    7  : problem 7
    8  : problem 8
    9  : problem 9
    10 : problem 10
%}
clc;
select = 8;

C_ox = 2.5e-15;
V_thn = .8;
k_n = 100e-6;
R_n = 20e3;
V_thp = .9;
k_p = 50e-6;
R_p = 40e3;
gnd = 0;


%------------------------------------------------------------------------------------
if (select == 1)
    vdd = 5;
    r1 = 100e3;
    r2 = r1;
    w_n = 6e-6;
    l_n = 0.6e-6;
    w_p = 6e-6;
    l_p = .6e-6;
    
    
    nmos_Voh = vdd;
    fprintf("\ntrivial, if NMOS gate gnd, ouputs: %0.2f\n", nmos_Voh);
    
    pmos_Vol = gnd;
    fprintf("trivial, if PMOS gate vdd, ouputs: %0.2f\n", pmos_Vol);
    
    %use square laws for Vsp...gate tied together == SATURATION
    syms vsp;
    nmos_vspL = (vdd - vsp) / r1;
    nmos_vspR = (k_n / 2) * (w_n / l_n) * (vsp - V_thn)^2;
    V_spn = double(solve( nmos_vspL == nmos_vspR, vsp));
    fprintf("\nV_spn possibilties:  %0.2f, %0.2f\n", V_spn(1), V_spn(2));
    fprintf("we must have V_gd <= V_thn, select V_spn = %0.2f\n", V_spn(2));
    
    pmos_vspL = vsp / r2;
    pmos_vspR = (k_p / 2) * (w_p / l_p) * (vdd - vsp - V_thp)^2;
    V_spp = double(solve(pmos_vspL == pmos_vspR, vsp));
    fprintf("\nV_spp possibilities:  %0.2f, %0.2f\n", V_spp(1), V_spp(2));
    fprintf("we must have V_dg <= abs(V_thp), select V_spp = %0.2f\n", V_spp(2));
    
    %nmos Vol, put 5V on gate, output goes low ... this is tridoe
    syms v_ol;
    nmos_volL = (vdd - v_ol) / r1;
    nmos_volR = k_n * (w_n / l_n) * (vdd - V_thn)*v_ol;
    nmos_vol = solve(nmos_volL == nmos_volR, v_ol);
    fprintf("nmos vol possibilities:  %0.4f\n", nmos_vol);
    
    syms voh;
    pmos_vohL = voh/r2;
    pmos_vohR = k_p * (w_p/l_p) * (vdd - V_thp) * (vdd-voh);
    pmos_voh = double(solve(pmos_vohL==pmos_vohR, voh))
    
    
end

%------------------------------------------------------------------------------------
if (select == 2)
    w_n = 6e-6;
    l_n = 0.6e-6;
    w_p = 6e-6;
    l_p = .6e-6;
    r1 = 100e3;
    r2 = r1;
    cap = 50e-15;
    
    nmos_tplh = 0.7 * r1 * cap;
    fprintf("\nNMOS t_PLH =  %0.2f  ns\n", nmos_tplh*10^9);
    nmos_tphl = 0.7 * (l_n / w_n) * R_n * cap;
    fprintf("NMOS t_PHL =  %0.2f  ns\n", nmos_tphl*10^9);
    
    pmos_tphl = 0.7 * r2 * cap;
    fprintf("\nPMOS t_PHL =  %0.2f  ns\n", pmos_tphl*10^9);
    pmos_tplh = 0.7 * (l_p / w_p) * R_p * cap;
    fprintf("PMOS t_PLH =  %0.2f  ns\n", pmos_tplh*10^9);
end

%------------------------------------------------------------------------------------
if (select == 3)
    dt = 100e-12;
    pd = 10e-9;
    vdd = 5;
    imx = 100e-3;
    
    %only pulls on an edge
    C_L = (imx * dt) / vdd;
    fprintf("\nusing i=c(dv/dt), CL = %0.2f pF\n", C_L*10^12);
    
    % capacitor doesn't disipate power
    I_avg = (C_L * vdd) / pd;
    P_avg = I_avg * vdd;
    fprintf("\navg pow: %0.3f W\n", P_avg);
    fprintf("PMOS dissipates half:  %0.2f mW\n", 1000*P_avg/2);
    fprintf("NMOS dissipates half:  %0.2f mW\n", 1000*P_avg/2);
    
    % power delivered to cap laod
    ecap = (1/2)*C_L*vdd^2;
    fprintf("\ncap gets energy:  %0.2f  pJ\n", ecap*10^12);
end

%------------------------------------------------------------------------------------
if (select == 4)
    wn1 = 30;
    ln1 = .6;
    wn2 = 90;
    ln2 = .6;
    c_1 = wn1 * ln1 * C_ox;
    c_2 = c_1;
    c_3 = (3/2) * wn2 * ln2 * C_ox;
    c_4 = c_3;
    ctot = c_1 + c_2 + c_3 + c_4;
    fprintf("\ntotal cap:  %0.2f fF\n", ctot*10^15);
    
    rnmos = (ln1/wn1)*R_n;
    rpmos = (ln1/wn1)*R_p;
    tplh = .7 * rpmos * ctot;
    tphl = .7 * rnmos * ctot;
    fprintf("t_PLH = %f  ps\n", tplh*10^12);
    fprintf("t_PHL = %f  ps\n", tphl*10^12);
    
    
end

%------------------------------------------------------------------------------------
if (select == 7)
    cap = 1e-12;
    capp = .25 * cap;
    v_in = 5;
    v_out = v_in * (4/5);
    fprintf("\nthe change is %0.2f V\n", v_out);
end

%------------------------------------------------------------------------------------
if (select == 8)
    cap = 1.25e-12;
    res = 20e3;
    
    syms ff;
    eqn = (2/1.25) * res * 1/sqrt(res^2 + (1/(2*pi*ff*cap))^2);
    freq = double(solve(eqn == 1, ff));
    fprintf("\nhis way:  %0.3f  MHz\n", freq(2)/(10^6));
    
    %eqq = sqrt(res^2 + (1/(2*pi*ff*cap))^2);
    %check = double(solve(eqq == 32e3, ff));
    
    eqq = res^2 + (1/(2*pi*ff*cap))^2;
    check = double(solve(eqq == (32e3)^2, ff));
    temp = (32e3)^2 - res^2;
    
    eqq = sqrt(1/((cap*2*pi)^2 * temp));
    ttt= pi
    tttt = 22/7
    
end

%------------------------------------------------------------------------------------
if (select == 9)
    
end

%------------------------------------------------------------------------------------
if (select == 10)
    
end
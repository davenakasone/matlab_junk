%{
    1  :  try out f_perm_ep() 
    2  :  try f_cmil2m(), f_m2cmil()
    3  :  try f_line_rp()
    4  :  try f_L_line_single()
    5  :  try f_L_line_three()
    6  :  try f_XL_3_meter(), f_XL_3_mile()
    7  :  try f_C_line_single(), f_C_line_three()
    8  :  try f_YC_line_three() , f_XC_line_three()
    9  :  pp9.1, find resistance in ohms/km
    10 :  try f_miles2km(), f_km2miles()
    11 : pp9.11, use f_line_mid_ABCD()
%}
close all;
clear all;
clc;
select = 12;


%-------------------------------------------------------------------------------------
if select == 1
    get_epp0 = f_perm_ep(1)
end


%-------------------------------------------------------------------------------------
if select == 2
    try_cmil = 1000; % diameter of wire in cmils
    to_meters = f_cmil2m(try_cmil)
    to_cmil = f_m2cmil(to_meters)
end


%-------------------------------------------------------------------------------------
if select == 3
    rr = 1;
    r_p = f_line_rp(rr)
end


%-------------------------------------------------------------------------------------
if select == 4
    D = 2;
    r = 3;
    L = f_L_line_single(D, f_line_rp(r)) % H/m
end


%-------------------------------------------------------------------------------------
if select == 5
    d_1 = 1;
    d_2 = 2;
    d_3 = 3;
    r = 4; 
    L = f_L_line_three(f_GMD(d_1, d_2, d_3), f_line_rp(r)) % H/m
end


%-------------------------------------------------------------------------------------
if select == 6
    d_1 = 50;
    d_2 = 60;
    d_3 = 70;
    rr = 1;
    XLmeter = f_XL_3_meter(f_GMD(d_1, d_2, d_3), f_line_rp(rr))  % ohms/meter
    XLmile = f_XL_3_mile(f_GMD(d_1, d_2, d_3), f_line_rp(rr))  % ohms/meter
end


%-------------------------------------------------------------------------------------
if select == 7
    rr = 0.1;
    D = 1;
    Cap1 = f_C_line_single(D, rr)  %  F/m
    Cap3 = f_C_line_three(f_GMD(1,2,3), rr) % F/m
end


%-------------------------------------------------------------------------------------
if select == 8
    d1=1;
    d2=2;
    d3=3;
    rr = 0.05;
    
    xx=f_YC_line_three(f_GMD(d1, d2, d3), rr)   % S/m
    xxx=f_XC_line_three(f_GMD(d1, d2, d3), rr)  % omh/m
    xxxx = 1/ xx
end


%-------------------------------------------------------------------------------------
if select == 9
    % given aluminum conductor with radius, find R_DC in ohms/km
    d = 3/100; % m
    r = d / 2; 
    rho = 2.83e-8; % from table
    R_DC_m = rho / (pi * r^2); % ohms per meter
    R_DC_km = 1000 * R_DC_m
end


%-------------------------------------------------------------------------------------
if select == 10
    km_z = 1000;
    m_z = f_km2mile(km_z)
    km_z = f_mile2km(m_z)
end


%-------------------------------------------------------------------------------------
if select == 11
    d_km = 100; %km
    % 3 phase
    V_line = 138e3; % V
    S_out = 200e6; % VA
    f = 60; % Hz
    w = 2 * pi * f; % rad/sec
    r_km = 0.103; % omhs / km
    x_km = 0.525; % ohms / km
    y_km = 3.3e-6; % S / km
    
    %a), per-phase series impedance and shund admittance of this line
    z_phase = (r_km + 1j*x_km) * d_km;
    f_mdri("Z_phase", z_phase, 1); 
    y_phase = 1j*y_km * d_km;
    f_mdri("Y_phase", y_phase, 1); 
    
    %b), should it be short, medium, or long?
    d_mi = f_km2mile(d_km); 
    fprintf("\nline is:  %0.3f  miles, use medium\n", d_mi);
    
    %c), find ABCD for network
    [A, B, C, D] = f_line_med_ABCD(z_phase, y_phase);
    f_mdri("A", A, 1);
    f_mdri("B", B, 1);
    f_mdri("C", C, 1);
    f_mdri("D", D, 1);
    
    %d), phasor digram when supplying rating @ pf = 0.9 lag
    
    %e), find Vs if line is supplying rating at S , pf = 0.9 lag
    V_phase = V_line / sqrt(3) % the rated line voltage reduces
    I_line = S_out / (sqrt(3) * V_line)
    V_r = V_phase % assume 0 deg
    
    %f), what is the voltage regulation
    %g), what is the efficency?
end

%-------------------------------------------------------------------------------------
if select == 12
    
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END>  
%{
    1  :  ex1, a circle of cast-steel, BH curve, and some turns
    2  :  ex2, cast iron and sheet steel on rectangle
    3  :  p1, core and air gap
    4  :  p2, circle, coil, air gap
    5  :  p3, force
    6  :  p4, force on spring  ???
    7  :  p5, block
    8  :  pp1.5, block
    9  :  pp1.6, 2 air gaps
    10 :  pp1.8, core w/ 3 legs
    11 :  pp1.10, wire in mag field, Lenz on CP
    12 :  pp1.12, 2 coils on block
    13 :  pp1.14, airgap on coil
    14 :  hw, assgn2, part 1
    15 :  hw, assgn2, part 2

    f_perm_mu()
%}
clear all;
clc;
select = 15;

%-------------------------------------------------------------------------------------
if select == 1
    A = 2e-3; % m^2, cross-section area of the ring
    l = 0.16; % m, the mean length of the ring
    N = 400; % 400 turns on the ring
    
    % what value of current I will make magnetic flux 0.4 mWb?
    phi = 4e-4; % given
    B = phi / A; % because phi = BA, solve for B
    H = 170; % because cast-steel curve says H = 170 @ B = 0.2
    I = H * l / N;
    fprintf("\n\tneed I = %0.1f mA\n", 1000*I);
    
    % determine mu of material under these conditions
    mu = B * l / (N * I);
    fprintf("\n\tpermiability is:  %0.7f\n", mu);
    mu_r = mu / f_perm_mu(1); % from  mu=mu_r * mu_0
    fprintf("\n\trelative perimiability:  % 0.1f\n", mu_r);
end

%-------------------------------------------------------------------------------------
if select == 2
    % 2 metals in core, want flux phi = 2.5e-4 Wb
    l_steel = ((4 + 4 + 4) * 2.54) / 100; % inches to cm to m
    l_iron = ((0.5 + 4 + 0.5) * 2.54) / 100; % inches to cm to m
    A = (2.54/100)^2; % in^2 to cm^2 to m^2
    phi = 3.5e-4; % given
    N = 50; % given
    B = phi / A; % because phi = BA, therefore B = 0.543 T
    H_steel = 70; % from BH curve for sheet steel
    H_iron = 1600; % from BH curve for cast iron
    I = (H_steel * l_steel + H_iron * l_iron) / N; % apply law
    fprintf("\n\tI = %0.2f A, by Ampere's Law\n", I);
end

%-------------------------------------------------------------------------------------
if select == 3
    flux = 2.4e-4; % Wb
    l_air = 0.003; % m
    l_ab = 0.05;
    l_ef = l_ab;
    l_af = 0.02;
    l_be = 0.02;
    l_bc = (l_be - l_air) / 2;
    l_de = l_bc;
    l_steel = l_ab + l_bc + l_de + l_ef;
    N = 100;
    A = 2e-4;
    B = flux / A;
    H_steel = 340; 
    H_air = B / f_perm_mu(1);
    Hl_steel = H_steel * l_steel
    Hl_air = H_air * l_air
    I = (Hl_steel + Hl_air) / N
    mu_r_steel = B / (H_steel * f_perm_mu(1))
end

%-------------------------------------------------------------------------------------
if select == 4
    phi = 2e-4; % flux you want
    N1 = 200;
    N2 = 40;
    I2 = 0.3;
    r = 0.3;
    circ = 2*pi*r;
    l_air = 0.002;
    l_steel = circ - l_air;
    A = 1.3e-4;
    
    B = phi / A; % matching sheet steel at B = 1.54 T
    H_steel = 300;
    H_air = B / f_perm_mu(1);
    Hl_air = H_air * l_air;
    Hl_steel = H_steel * l_steel;
    
    mur = B / (H_steel * f_perm_mu(1));
    phi2 = N2*I2*( (mur*f_perm_mu(1)*A)/l_steel + (f_perm_mu(1)*A)/l_air)
    phi1 = phi - phi2
    I1 = (phi1/N1) * ( l_steel/(mur*f_perm_mu(1)*A) + l_air/(f_perm_mu(1)*A) )
end

%-------------------------------------------------------------------------------------
if select == 5
    x = 4e-2;
    N = 80;
    I = 0.9;
    d_phi = (8e-4 - 0.5e-4) / (x/2); % rate of change, flux vs dist
    f = (1/2) * N * I * d_phi
end

%-------------------------------------------------------------------------------------
if select == 6
    D = 0.01;
    A = 2 * pi * (D/2);
    N = 200;
    l_air = 0.2/100;
    phi = 0.2e-4;
    
    I = (phi/N) * (l_air/(f_perm_mu(1)*A))
    B = N * I * f_perm_mu(1) / l_air
    f = (1/2) * (B^2 * A / f_perm_mu(1))
end

%-------------------------------------------------------------------------------------
if select == 7
    N = 50;
    phi = 1.5e-4;
    l = 0.2;
    A = 6e-4;
    phi2 = 1.5e-4;
    phi1 = phi;
    phiT = phi1 + phi2;
    
    B = phi/A;
    H_steel = 40;
    %I = (1/N) * (
end

%-------------------------------------------------------------------------------------
if select == 8
    phi = 0.003;
    depth = 5/100;
    N = 500;
    l_left = ( (15/2) + 15 + (15/2) ) / 100; % coils are here
    A_left = (10/100) * depth;
    l_tb = 2 * ( (10/2) + 20 + (5/2) ) / 100; % the top and bottom taken as 1
    A_tb = (15/100) * depth;
    l_right = l_left;
    A_right = (5/100) * depth;
    mu_r = 1000;
    mu = mu_r * f_perm_mu(1);
    
    rel_tb = l_tb / (mu * A_tb);
    rel_left = l_left / (mu * A_left);
    rel_right = l_right / (mu * A_right);
    rel_tot = rel_tb + rel_left + rel_right;
    I = phi * rel_tot / N;
    fprintf("\n\tneed %0.2f A \n", I);
    
    B_tb = phi / A_tb;
    fprintf("\n\tB on top/bottom = %0.2f T\n", B_tb);
    B_left = phi / A_left;
    fprintf("\n\tB on left = %0.2f T\n", B_left);
    B_right = phi / A_right;
    fprintf("\n\tB on right = %0.02f T\n", B_right);
end

%-------------------------------------------------------------------------------------
if select == 9
    mur = 2000;
    mu = mur * f_perm_mu(1);
    depth = 7/100;
    A = depth^2;
    l_agl = 0.05/100;
    l_agr = 0.07/100;
    N = 300;
    I = 1;
    l_center = (7/2 + 30 + 7/2) / 100;
    l_left = ((30 + 7 + 30 + 7 + 30 + 7)/100) - l_agl;
    l_right = ((30 + 7 + 30 + 7 + 30 + 7)/100) - l_agr;
    
    r_agl = l_agl / (f_perm_mu(1)*A*1.05);
    r_agr = l_agr / (f_perm_mu(1)*A*1.05);
    r_center = l_center / (mu * A);
    r_left = l_left / (mu*A);
    r_right = l_right / (mu*A);
    r_tot = r_center + f_para(r_agr + r_right, r_agl + r_left);
    
    phi = N * I / r_tot
    temp = (r_agl + r_left) / (r_agl + r_left + r_agr + r_right);
    phi_left = temp * phi
    B_left = phi_left / (A*1.05)
    temp = (r_agr + r_right) / (r_agl + r_left + r_agr + r_right);
    phi_right = temp * phi
    B_right = phi_right / (A*1.05)
end

%-------------------------------------------------------------------------------------
if select == 10
    depth = 5/100;
    N = 200;
    mu = f_perm_mu(1500);
    mu0 = f_perm_mu(1);
    A_lr = (9/100) * depth;
    A_center = (15/100) * depth;
    l_ag = 0.04/100;
    I = 2;
    l_left = (37 + 34 + 37) / 100;
    l_center = ((25 + 9) / 100) - l_ag;
    l_right = l_left;
    
    r_ag = l_ag / (mu0*A_center*1.04);
    r_left = l_left / (mu*A_lr);
    r_right = r_left;
    r_center = l_center / (mu*A_center);
    r_tot = r_left + f_para(r_center + r_ag, r_right);
    phi = (N*I)/r_tot;
    fprintf("\n\ttotal flux is from left leg:  %0.5f Wb\n", phi);
    temp = r_right / (r_right + r_center + r_ag);
    phi_center = temp * phi;
    fprintf("\n\tflux in center:  %0.5f\n", phi_center);
    temp = (r_center + r_ag) / (r_right + r_center + r_ag);
    phi_right = temp * phi;
    fprintf("\n\tflux in right:  %0.5f\n", phi_right);
    check = phi - phi_center - phi_right
end

%-------------------------------------------------------------------------------------
if select == 11
    v = 6; % m/s
    l = 0.75; % m
    tht = pi/4; % 45 degrees
    B = 0.2; % T into page
    vel = [-v, 0, 0];
    BB = [0, 0, -B];
    v_cp_B = cross(vel, BB); % going down Y axis
    e_ind = v*B*l*cos(tht); % 0.6364 volts down....
end

%-------------------------------------------------------------------------------------
if select == 12
    depth = 15/100;
    A = depth^2;
    N1 = 600;
    N2 = 200;
    l = (50 + 15 + 50 + 15 + 50 + 15 + 50 + 15)/100;
    I1 = 0.5;
    I2 = 1;
    F_tot = N1*I1 + N2*I2;
    H = F_tot / l; % match this H = 192.3
    B = 0.16; % from given curve
    phi = B * A
    mu = (phi*l)/(A*F_tot);
    mur = mu/f_perm_mu(1);
end

%-------------------------------------------------------------------------------------
if select == 13
    depth = 5/100;
    A_ag = depth^2;
    l_ag = 0.07/100;
    N = 500;
    B = 0.5;
    l = (2*(40 + 37.5)-l_ag)/100;
    
    phi = B*A_ag*1.05; % with excess 5%
    B_right = phi/A_ag;
    B_other = phi/(depth*(10/100)); % other wider parts
    H_ag = B/f_perm_mu(1);
    H_right = 410;
    H_other = 240;
    F_tot = H_ag*l_ag + H_right*0.40 + H_other*3*0.40;
    i = F_tot/N; 
end

%-------------------------------------------------------------------------------------
if select == 14
    mu = 5e-3;
    mur = 4000;
    mu0 = f_perm_mu(1);
    A = 1e-2;
    V = 120; % RMS
    f = 60;
    w = 2*pi*f;
    N = 400;
    l_left = 3*0.24/100;
    l_right = l_left;
    l_center = 0.24/100;
    
    rel_left = l_left / (mu * A);
    rel_right = rel_left;
    rel_center = l_center / (mu * A);
    rel_tot = rel_center + f_para(rel_left, rel_right);
    fprintf("\n\treluctance seen from source:  %0.1f  At/Wb\n", rel_tot);
    
    phi = V/(w*N);
    fprintf("\n\tPhi rms = %0.3f mWb\n", phi*1000);
    
    L = N^2 / rel_tot;
    fprintf("\n\tL=  %0.3f  H\n", L);
    
    I = phi*rel_tot/N;
    fprintf("\n\tI=  %0.3f mA\n", I*1000);
    
    phi_right = (rel_left/(rel_left + rel_right)) * phi;
    fprintf("\n\tphi_right = %0.3f mWb\n", 1000*phi_right);
    E_right = phi_right * 100 * w;
    fprintf("\n\tthe new coil potentail= %0.2f Vrms\n", E_right);
end

%-------------------------------------------------------------------------------------
if select == 15
    mu = 5e-3;
    mur = 4000;
    mu0 = f_perm_mu(1);
    A = 1e-2;
    V = 120; % RMS
    f = 60;
    w = 2*pi*f;
    N = 400;
    l_left = 3*0.24/100;
    l_right = l_left;
    l_center = 0.24/100;
    
    rel_left = l_left / (mu * A);
    rel_right = rel_left;
    rel_center = l_center / (mu * A);
    rel_tot = rel_left + f_para(rel_center, rel_right);
    fprintf("\n\treluctance seen from source:  %0.1f  At/Wb\n", rel_tot);
    
    phi = V/(w*N);
    fprintf("\n\tPhi rms = %0.3f mWb\n", phi*1000);
    
    L = N^2 / rel_tot;
    fprintf("\n\tL=  %0.3f  H\n", L);
    
    I = phi*rel_tot/N;
    fprintf("\n\tI=  %0.3f mA\n", I*1000);
    
    phi_right = (rel_center/(rel_center + rel_right)) * phi;
    fprintf("\n\tphi_right = %0.3f mWb\n", 1000*phi_right);
    E_right = phi_right * 100 * w;
    fprintf("\n\tthe new coil potentail= %0.2f Vrms\n", E_right);
end
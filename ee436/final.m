%{
    ee436 final exam
    1 : problem 1
    2 : problem 2
%}
clc;
close;
clear;
select = 1;


if select == 1
    s = (1/sqrt(2)) .* ...
        [
            0 , 0                         , 1                         , 1;
            0 , 0                         , 1 * exp(1j * deg2rad(180)), 1;
            1 , 1 * exp(1j * deg2rad(180)), 0                         , 0;
            1 , 1                         , 0                         , 0
        ];
    display(s);
    fprintf("\n [s] - transpose([s])=\n\n")
    display(s - transpose(s));
    fprintf("\nconj([s]) * transpose([s]) \n\n");
    display(conj(s) * transpose(s));
    fprintf("\n[s][s.T]=\n");
    display(s * transpose(s));

p_in = 1; 
g1 = 0.5;
g2 = 0.6;
g4 = 0.8;
gden = (1 - g1*g4) + (1 - g2*g4);
v1n3 = sqrt(2) * (1 - g2*g4) / gden;
v2n3 = -1 * sqrt(2) * (1 - g1*g4) / gden;
v3n3 = ((1 - g2*g4)*g1 + (1 - g1*g4)*g2) / gden;
v4n3 = (g1 - g2) / gden;

v3p = sqrt((2 * p_in) / (1 - v3n3^2));
v3n = v3p * v3n3;
g3 = v3n / v3p;
v1n = v3p * v1n3;
v2n = v3p * v2n3;
v3n = v3p * v3n3;
v4n = v3p * v4n3;
v1p = g1 * v1n;
v2p = g2 * v2n;
v4p = g4 * v4n;

p1avg = (1/2) * (abs(v1n)^2 - abs(v1p)^2);
p1l = p1avg * (1-g1^2);

p2avg = (1/2) * (abs(v2n)^2 - abs(v2p)^2);
p2l = p2avg * (1-g2^2);

p4avg = (1/2) * (abs(v4n)^2 - abs(v4p)^2);
p4l = p4avg * (1-g4^2);

p3l = p1l + p2l + p4l;

fprintf("\n~V1n  =  %9.3f  * ~V3^+  =  %9.3f,  ~V1p=  %9.3f  ,  G1=  %0.2f (%0.2f)\n", v1n3, v1n, v1p, v1p/v1n, g1);
fprintf("~V2n  =  %9.3f  * ~V3^+  =  %9.3f,  ~V2p=  %9.3f  ,  G2=  %0.2f (%0.2f)\n"  , v2n3, v2n, v2p, v2p/v2n, g2);
fprintf("~V3n  =  %9.3f  * ~V3^+  =  %9.3f,  ~V3p=  %9.3f  ,  G3=  %0.2f (%0.2f)\n"  , v3n3, v3n, v3p, v3n/v3p, g3);
fprintf("~V4n  =  %9.3f  * ~V3^+  =  %9.3f,  ~V4p=  %9.3f  ,  G4=  %0.2f (%0.2f)\n"  , v4n3, v4n, v4p, v4p/v4n, g4);

fprintf("\np_in = %d W = (1/2) * ( |V3p|^2 - |V3n|^2 )  = %0.2f W\n",...
    p_in, 0.5 * (abs(v3p)^2 - abs(v3n)^2));

fprintf("\nP_1L=  %0.4f W, power supplied to port#1\n", p1l);
fprintf("P_2L=  %0.4f W, power supplied to port#2\n", p2l);
fprintf("P_4L=  %0.4f W, power supplied to port#4\n", p4l);
fprintf("P_3D=  %0.4f W, |P_1L + P_2L + P3L|, power absorbed by port#3\n\n", p3l);


end


%%%%~~~~

    
if select == 2
    % given:
    Z_L = 30;
    Z_S = 20;
    Z_0 = 50;
    s = [...
        0.45 * exp(1j * deg2rad(150)), 0.01 * exp(1j * deg2rad(-10));
        2.05 * exp(1j * deg2rad(10)) , 0.4  * exp(1j * deg2rad(-150))...
        ];
    
    % find the Gammas:
    gam_S = (Z_S - Z_0) / (Z_S + Z_0);
    gam_L = (Z_L - Z_0) / (Z_L + Z_0);
    gam_in = s(1,1) + (s(1,2) * s(2,1) * gam_L) / (1 - s(2,2) * gam_L);

    % calculate G_T :
    gt_num_a = abs(s(2,1))^2;
    gt_num_b = 1 - abs(gam_S)^2;
    gt_num_c = 1 - abs(gam_L)^2;
    gt_num = gt_num_a * gt_num_b * gt_num_c;
    gt_den_a = abs(1 - gam_S * gam_in)^2;
    gt_den_b = abs(1 - s(2,2) * gam_L)^2;
    gt_den = gt_den_a * gt_den_b;
    GT = gt_num / gt_den;
    GTdb = 10 * log10(GT);

    display(s);
    fprintf("\n[s] - transpose([s]) = \n\n");
    display(s - transpose(s));
    fprintf("\nconj([s]) * transpose([s]) = \n\n");
    display(conj(s) * transpose(s));

    fprintf("\nGamma_S = %0.3f\n", gam_S);
    fprintf("Gamma_L = %0.3f\n", gam_L);
    p_cplx("Gamma_in", gam_in);
    fprintf("\nGT = %0.3f     ...%0.3f dB\n\n", GT, GTdb);
end


%%%%~~~~


if select == 99
    fprintf("\n\n\t\tdone\n\n");
end


%%%%~~~~END>  final.m





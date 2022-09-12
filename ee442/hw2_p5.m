% hw2, p5, work problem 3-17, with R=1k
clc;

f = 60;
R = 1e3;
C = 1000e-6;
RC = R*C;
w = 2 * pi * f;
V_m = 200;

%a)
T = 1/f;
ratio_rct = RC / T;
fprintf("a) ratio RC/T = %0.1f\n", ratio_rct);
fprintf("a) this is big, so output voltage will decay slow, relative to input, when diode is off\n");

%b) exact ripple
theta = atan(-1 * w * R * C);
syms alp;
eqn_a = sin(alp);
eqn_b = sin(theta) * exp( (-1 * (2*pi + alp - theta) ) / (w * R * C));
alpha = vpasolve(eqn_a == eqn_b, alp);
check = double(subs(eqn_a - eqn_b, alp, alpha));
fprintf("\nchecked:  %0.3f, should be 0\n", check);
vripple_exact = V_m * (1 - sin(-1*alpha));
fprintf("b) exact peak-peak ripple voltage:  %0.4f V\n", vripple_exact);

%c) approx ripple
vripple_approx = V_m / (f*R*C);
fprintf("c) approximate peak-peak ripple:  %0.4f V\n", vripple_approx);
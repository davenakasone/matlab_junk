% ex9.6, design a lag-lead compensator, page 874
format compact;
clear;
close all;
clc;

% step 1, the uncompensated system
o_s = 0.2; % need 20% overshoot
zeta = f_os_2_zeta(o_s);
fprintf("\ngiven the OS, zeta must be:  %0.3f\n", zeta);
theta = pi - acos(zeta);
fprintf("\nfinding zeta, theta is:  %0.3f rad  ,  %0.3f deg\n",...
    theta, rad2deg(theta));

syms s;
syms K;
Gs_num = 1;
gsn = 1;
Gs_den = s * (s + 6) * (s + 10);
gsd = sym2poly(Gs_den);
Gs = Gs_num / Gs_den;
%pretty(expand(Gs));
tf_raw = tf(gsn, gsd)
%f_rlocus(tf_raw);
K_raw = 192.1;
s_raw = -1.794 + 1j*3.501;
fprintf("\nbook wants to use:  K = %0.1f  ,  s = %0.3f +/- j %0.3f\n",...
    K_raw, real(s_raw), imag(s_raw));
fprintf("\nnote the angle:  %0.3f rad  ,  %0.3f deg == theta\n",...
    angle(s_raw), rad2deg(angle(s_raw)));

% step 2, design the lead and handle settling time
new_wnz = 2 * real(s_raw); 
new_wd = new_wnz * tan(theta);
fprintf("\nthe new [ zeta * w_n] = %0.3f  ,  distance on real axis\n", new_wnz);
fprintf("\nthe new [ w_d ] = %0.3f  ,  distance on imaginary axis\n", new_wd);

s_lead = new_wnz + 1j * new_wd;
v1 = f_make_vec(0, s_lead);
v2 = f_make_vec(-10, s_lead);
theta_pc_lead = pi - (angle(v1) + angle(v2));
fprintf("\nthe angle of the pole must be:  %0.3f rad  ,  %0.3f deg\n",...
    theta_pc_lead, rad2deg(theta_pc_lead));

syms pc_l;
eqn_pc_lead = imag(s_lead) / (pc_l - abs(real(s_lead)));
pc_lead = double(solve(eqn_pc_lead==tan(theta_pc_lead), pc_l));
fprintf("\nthe compensated pole is located at pc_lead=  %0.3f\n", pc_lead);

%steps 3 and 4
KGs = KK * Gs;
G_lead = (s + 6) / (s + pc_lead);
new_Gs = G_lead * KGs;
[new_Gs_num, new_Gs_den] = numden(new_Gs);
new_gsn = sym2poly(new_Gs_num);
new_gsd = sym2poly(new_Gs_den);
tf_lead = tf(new_gsn, new_gsd)
f_rlocus(tf_lead);
T_raw = KGs / (1 + KGs);
[T_raw_num, T_raw_den] = numden(T_raw);
T_tf_raw = tf(sym2poly(T_raw_num), sym2poly(T_raw_den));
T_lead = new_Gs / (1 + new_Gs);
[T_lead_num, T_lead_den] = numden(T_lead);
T_tf_lead = tf(sym2poly(T_lead_num), sym2poly(T_lead_den));
%
figure('Position', [20, 20, 800, 800]);
hold on;
grid on;
step(T_tf_raw, 'r-', T_tf_lead, 'b--');
set(findall(gcf,'type','line'),'linewidth',2);
legend('raw','lead compensated','Location','best', 'fontsize', 14);
hold off;
%

% step 5
kv_raw = subs(s*KGs, s, 0);
e_raw = 1 / kv_raw;
kv_lead = subs(s*new_Gs, s, 0);
e_lead = 1 / kv_lead;
fprintf("\nerror, no compensation=    %0.3f  ,  Kv=  %0.3f\n", e_raw, kv_raw);
fprintf("error, lead compensation=  %0.3f  ,  Kv=  %0.3f\n", e_lead, kv_lead);

%%%%%%%%~~~~~~~~END>
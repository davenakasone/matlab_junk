%{
    1  :  data table
%}
clc;
close all;
clear;

n_min = 3;
n_max = 9;
data_sets = n_max - n_min + 1;
L_m = 1.3631; % FILL
L_cm = round(L_m * 100, 5, "significant");
rho_s = 6.83e-4; % FILL    kg/m
f = 120; % FILL
fprintf("\nrho_string:  %d  kg/m\n", rho_s);
fprintf("L = %0.4f m  ,  %0.2f cm\n\n", L_m, L_cm);


mass_g = [785, 445, 285, 200, 145, 110, 90]; % FILL
mass_kg = zeros(1, data_sets);

fprintf("n   |  mass(g)  |   mass(kg)\n");
fprintf("-----------------------------\n");
for ii = 1:1:data_sets
    mass_kg(1, ii) = round(mass_g(1, ii) / 1000, 3, "significant");
    fprintf("%d   |  %3d      |  %0.3f\n", ii+2, mass_g(1, ii), mass_kg(1, ii));
end

tension = 9.8 .* mass_kg;
lambda_cm = [90.873, 68.155, 54.524, 45.437, 38.946, 34.078, 30.291]; % FILL
lambda_m = lambda_cm ./ 100;
tension_sqrt = sqrt(tension);
fprintf("\nT=Mg   |  lam_cm  |  lam_m    |  sqrt(T)\n");
fprintf("------------------------------------------\n");
for ii = 1:1:data_sets
    fprintf("%0.2f   |  %0.3f  |  %0.5f  |  %0.3f\n",...
        tension(1, ii), lambda_cm(1, ii), lambda_m(1, ii), tension_sqrt(1, ii));
end

m_a  = data_sets .* dot(tension_sqrt, lambda_m);
m_b = sum(tension_sqrt) * sum(lambda_m);
m_c = data_sets .* sum(tension);
m_d = sum(tension_sqrt)^2;
mm = (m_a - m_b) / (m_c - m_d);

b_a = sum(lambda_m) * sum(tension);
b_b = dot(tension_sqrt, lambda_m) * sum(tension_sqrt);
b_c = data_sets .* sum(tension);
b_d = sum(tension_sqrt)^2;
bb = (b_a - b_b) / (b_c - b_d);

r_a = data_sets * dot(tension_sqrt, lambda_m);
r_b = sum(tension_sqrt) * sum(lambda_m);
r_c = data_sets * sum(tension) - sum(tension_sqrt)^2;
r_d = data_sets * sum(lambda_m.^2) - sum(lambda_m)^2;
rr = (r_a - r_b) / (sqrt(r_c) * sqrt(r_d));

fprintf("\nY= mx + b:::   Y = [%0.4f]x + [%0.4f]  ,  r=  %0.5f\n\n", mm, bb, rr);
f_exp = 1/(mm*sqrt(rho_s));
f_error = 100 * abs(f_exp-f) / f;
fprintf("#1 f_exp=  %0.1f  Hz,  error:  %0.1f %%  ...good accuracy\n", f_exp, f_error);

xx = linspace(min(tension_sqrt,[], "all"), max(tension_sqrt, [], "all"), 100);
Y =  mm .* xx + bb;
figure()
hold on;
grid on;
axis padded;
title("\lambda vs sqrt(Tension)", FontSize=20);
xlabel("sqrt(Tension) in N", FontSize=15);
ylabel("\lambda in m", FontSize=15);
plot(tension_sqrt, lambda_m, "rx", MarkerSize=10, LineWidth=2);
plot(xx, Y, "b-", LineWidth=1);
hold off;

fprintf("#2\n")
fexp_array = (1 ./ lambda_m) .* sqrt(tension ./ rho_s);
fexp_avg = sum(fexp_array) / data_sets;
std_dev = 0;
for ii = 1:1:data_sets
    fprintf("n=%d  ,  f=  %0.1f\n", ii+2, fexp_array(ii));
    std_dev = std_dev + (fexp_array(ii) - fexp_avg)^2;
end
std_dev = sqrt((1/(data_sets-1)) * std_dev);
std_err = std_dev / sqrt(data_sets);
fprintf("f_exp_avg:  %0.1f\n", fexp_avg);
fprintf("std err: %0.2f\n", std_err);

fprintf("\n#3\n");
lam2 = 2 * L_m / 2;
lam8 = 2 * L_m / 8;
v2 = round(lam2 * f, 3, "significant");
v8 = round(lam8 * f, 3, "significant");
fprintf("n=2, v=  %0.1f\n", v2);
fprintf("n=8, v=  %0.1f\n", v8);

fprintf("\n#4\n");
lam1 = 2 * L_m / 1;
t1 = rho_s * (f * lam1)^2;
fprintf("T=  %0.1f N\n", t1);



%%%%%%%%%~~~~~~~~END>  lab21.m

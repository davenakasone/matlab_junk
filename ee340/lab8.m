clc;
clear all;
close all;

rpm =     [1800   , 1775  , 1750  , 1725  , 1700  , 1675  , 1650  , 1625 , 1600  , 1575  , 1550  ];
torque =  [0.007  , 0.355 , 0.687 , 0.985 , 1.253 , 1.493 , 1.701 , 1.893, 2.058 , 2.203 , 2.333 ];
P_mech =  [132.0  , 659.3 , 126.0 , 177.9 , 222.9 , 261.8 , 293.9 , 322.1, 344.9 , 363.4 , 378.7 ];
current = [0.774  , 0.788 , 0.85  , 0.944 , 1.057 , 1.182 , 1.303 , 1.431, 1.555 , 1.673 , 1.791 ];
P =  3 .* [12.76  , 34.76 , 56.91 , 77.66 , 97.54 , 116.7 , 133.8 , 150.7, 166.8 , 181.3 , 195.5 ];
Q =  3 .* [95.58  , 91.7  , 89.19 , 88.08 , 88.18 , 89.32 , 91.25 , 94.01, 97.34 , 101.1 , 105.4 ];
theta =      [-82.39 , -69.33, -57.57, -48.65, -42.18, -37.54, -34.35, -32  , -30.37, -29.23, -28.43];

temp = size(rpm);
data_length = temp(2);
eff = zeros(1, data_length);
p_f = zeros(1, data_length);
for ii = 1:1:data_length
    eff(1,ii) = P_mech(1,ii) / P(1,ii);
    p_f(1, ii) = cos(deg2rad(theta(1,ii)));
end
eff(1,1) = 0;
eff(1,2) = P(1,2)/P_mech(1,2);
fprintf("\nSpeed (RPM) | Torque (Nm) | P_mech (W) | Current (A)  |  Px3 (W) | Qx3 (VAR) | pf    | eta\n");
fprintf("---------------------------------------------------------------------------------------------\n");
for ii = 1:1:data_length
    fprintf("%7.2f     | ", rpm(1,ii));
    fprintf("%7.2f     | ", torque(1,ii));
    fprintf("%7.2f    | ", P_mech(1,ii));
    fprintf("%7.2f      | ", current(1,ii));
    fprintf("%7.2f  | ", P(1,ii));
    fprintf("%7.2f   | ", Q(1,ii));
    fprintf("%5.3f | ", p_f(1,ii));
    fprintf("%5.3f\n", eff(1,ii));
end

position_ = [20, 20, 400, 400];
figure(name=string(1), Position=position_)
hold on;
grid on;
title("tourqe vs RPM", FontSize=20, FontWeight="bold");
xlabel("speed (RPM)", FontSize=14);
ylabel("tourque (Nm)", FontSize=14);
plot(rpm, torque, "r-", LineWidth=3);
hold off;

figure(name=string(2), Position=position_)
hold on;
grid on;
title("efficency vs RPM", FontSize=20, FontWeight="bold");
xlabel("speed (RPM)", FontSize=14);
ylabel("eta", FontSize=14);
plot(rpm, eff, "b-", LineWidth=3);
hold off;

figure(name=string(3), Position=position_)
hold on;
grid on;
title("current vs RPM", FontSize=20, FontWeight="bold");
xlabel("speed (RPM)", FontSize=14);
ylabel("current (A)", FontSize=14);
plot(rpm, current, "g-", LineWidth=3);
hold off;

figure(name=string(4), Position=position_)
hold on;
grid on;
title("power factor vs RPM", FontSize=20, FontWeight="bold");
xlabel("speed (RPM)", FontSize=14);
ylabel("power factor", FontSize=14);
plot(rpm, p_f, "c-", LineWidth=3);
hold off;

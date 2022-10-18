%{
    lab 21
%}
clc;
close all;
clear;

d=5;
Vm = 20;
T_room = 26.5 + 273;
R = 8.31;

sz1 = 16;
pos1 = [5  ,6  ,7  ,8  ,9  ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20];
kpa1 = [352,300,262,231,207,188,170,158,146,133,124,120,111,105,100,95];
vol1 = zeros(1, sz1);
n1 = zeros(1,sz1);

sz2 = 18;
pos2 = [3  , 4 , 5  ,6  ,7  ,8  ,9  , 10, 11,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20];
kpa2 = [285,225, 183,154,133,120,107, 99, 89, 81,77,72,69,65,61,57,54,51];
vol2 = zeros(1, sz2);
n2 = zeros(1,sz2);

sz3 = 20;
pos3 = [1  , 2  , 3  , 4  , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
kpa3 = [338, 209, 152, 119, 98, 85, 74, 68, 61, 55, 50, 47, 44, 41, 38, 37, 35, 34, 32, 30];
vol3 = zeros(1, sz3);
n2 = zeros(1,sz3);

for ii = 1:1:sz1
    vol1(1,ii) = pos1(1,ii) * pi * (d^2/4) + Vm;
    n1(1,ii) = kpa1(1,ii) * vol1(1,ii) / (R * T_room);
    %fprintf("volume(%2d cm):  %9.1f\n", pos1(1,ii), vol1(1,ii));
end
for ii = 1:1:sz2
    vol2(1,ii) = pos2(1,ii) * pi * (d^2/4) + Vm;
    n2(1,ii) = kpa2(1,ii) * vol2(1,ii) / (R * T_room);
    %fprintf("volume(%2d cm):  %9.1f\n", pos2(1,ii), vol2(1,ii));

end
for ii = 1:1:sz3
    vol3(1,ii) = pos3(1,ii) * pi * (d^2/4) + Vm;
    n3(1,ii) = kpa3(1,ii) * vol3(1,ii) / (R * T_room);
    %fprintf("volume(%2d cm):  %9.1f\n", pos3(1,ii), vol3(1,ii));
end

fprintf("avg1, n=  %0.3f mMol\n", sum(n1)/sz1);
fprintf("avg2, n=  %0.3f mMol\n", sum(n2)/sz2);
fprintf("avg3, n=  %0.3f mMol\n", sum(n1)/sz3);

figure()
hold on;
grid on;
axis padded;
xlabel("V in cm^3", FontSize=15);
ylabel("p in kPa", FontSize=15);
plot(vol1, kpa1, 'rx' , linewidth=2);
plot(vol2, kpa2, 'bx' , linewidth=2);
plot(vol3, kpa3, 'gx' , linewidth=2);
hold off;



%%%%%%%%%~~~~~~~~END>  lab25.m

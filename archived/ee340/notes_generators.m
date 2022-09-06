%{
    synchronous generators
    field == rotor
    armeture == stator

    1  :  part1, slide 25, example 1  
%}
close all;
clear all;
clc;
select = 1;


%-------------------------------------------------------------------------------------
if select == 1
    % Y connected
    S_rate = 200e3;
    V_rate = 480;
    f = 50;
    If_rate = 5;
    VT_oc = 540;
    IL_sc = 300;
    R_armature = 0.2;
    
    % find model at rated conditions..reactance, resitance
    E_A = VT_oc / sqrt(3);
    fprintf("stator voltage is output Vf, open...  =  %0.3f V\n", E_A);
    
    % find the reactance of the armeture
    % Xs = sqrt(zs^2 - ra^2)
    X_s = sqrt( (E_A/IL_sc)^2 - R_armature^2);
    fprintf("stator reactance:  %0.3f ohm\n", X_s);
end


%-------------------------------------------------------------------------------------
if select == 99
 
end
%%%%%%%%~~~~~~~~~END> 

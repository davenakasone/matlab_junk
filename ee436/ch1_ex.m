%{
    chapter 1, EM

    1   :  ex1.1, plane wave, lossless
    2   :  ex1.

%}
clc;
clear;
select = 1;


if select == 1
    % plane wave, lossles dielectric medium
    % E_x = E_0 cos(wt - bz)
    % f = 5 GHz, wavelength = 3 cm
    f = 5e9;
    lamda = 3e-2;
    w = 2 * pi * f;

    k = 2 * pi / lamda; 
    fprintf("propogation constant  :  %0.3f  1/meter\n", k);

    v_p = w/k;
    fprintf("phase velocity        :  %0.3f  m/s   ...slower than speed of light\n", v_p);

    ep_r = (konst.c/v_p)^2;
    fprintf("realtive permittivity :  %0.3f  F/m\n", ep_r);

    eta_wave = konst.eta0/sqrt(ep_r);
    fprintf("wave impedance        :  %0.3f  ohms\n", eta_wave);
end

%%%%%%%%~~~~~~~~END>  ch1_ex.m

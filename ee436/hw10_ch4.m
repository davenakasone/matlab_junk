%{
    hw10, ch4
%}
clc;
close;
clear;
select = 1;


if select == 1
    s = [...
            0.178*exp(1j*deg2rad(90)), 0.6*exp(1j*deg2rad(45)) , 0.4*exp(1j*deg2rad(45)) , 0;
            0.6*exp(1j*deg2rad(45))  , 0                       , 0                       , 0.3*exp(1j*deg2rad(-45));
            0.4*exp(1j*deg2rad(45))  , 0                       , 0                       , 0.5*exp(1j*deg2rad(-45));
            0                        , 0.3*exp(1j*deg2rad(-45)), 0.5*exp(1j*deg2rad(-45)), 0
        ];
    U = s * conj(transpose(s));
    fprintf("\n|s| = \n\n");
    display(abs(s));
    display(U);

    
    display(s);
    fprintf("s.T =\n\n");
    display(transpose(s));
    fprintf("\n [s] - [s].T = \n\n");
    display(s-transpose(s));

    g = s(1,1) - s(1,3)*s(3,1);
    p_cplx("G1 ", g);
end
    


%%%%~~~~


if select == 99
    fprintf("\n\n\t\tdone\n\n");
end


%%%%~~~~END>  hw10_ch4.m





%{
    shunt admitance per phase @ 60 Hz
    S/m   =   2 pi f C
%}
function Y_C = f_YC_line_three(GMD, r_r)
    Y_C = (2.1 / 10^8) / log(GMD/r_r);
end
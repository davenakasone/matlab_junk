%{
    the inductance per phase in a 3-phase transmission line
    L in H/M
    GMD = (d1 d2 d3) ^ (1/3)
    GMR = 0.7788 r for solid
%}
function inductance = f_L_line_three(GMD, GMR)
    inductance = (2 / 10^7) * log(GMD/GMR);
end
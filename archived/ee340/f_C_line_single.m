%{
    capacitance in single phase
    C  in F/m
%}
function capacitance = f_C_line_single(D_D, r_r)
    temp_a = 2 * pi * f_perm_ep(1);
    temp_b = log(D_D/r_r);
    capacitance = temp_a / temp_b;
end
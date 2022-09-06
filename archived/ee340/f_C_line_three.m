%{
    capacitance per phase in 3-wire
    C  in F/m
%}
function capacitance = f_C_line_three(GMD, r_r)
    temp_a = 2 * pi * f_perm_ep(1);
    temp_b = log(GMD/r_r);
    capacitance = temp_a / temp_b;
end
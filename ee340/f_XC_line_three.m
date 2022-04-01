%{
    shunt capacitance ... reactance per phase
    ohm/m
%}

function reactance = f_XC_line_three(GMD, r_r)
    reactance = (47.7 * 10^6) / log(GMD/r_r);
end
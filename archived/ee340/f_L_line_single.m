%{
    for single phase transmission line
    solves for inductance L    in    H/m
%} 
function inductance = f_L_line_single(distance, radius_p)
    temp = log(distance/radius_p);
    inductance = (4 / 10^7) * temp;
end
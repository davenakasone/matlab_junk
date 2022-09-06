%{
    get the 2-port parameters, given total admittance and impedance
%}
function [AA, BB, CC, DD] = f_line_med_ABCD(z_total, y_total)
    AA = (z_total * y_total / 2) + 1;
    BB = z_total;
    CC = y_total * ((z_total * y_total / 4) + 1);
    DD = AA;
end
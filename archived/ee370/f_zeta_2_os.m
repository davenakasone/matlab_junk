function os_out = f_zeta_2_os(zeta_in)
    temp_a = -1 * zeta_in * pi;
    temp_b = sqrt(1 - zeta_in^2);
    os_out = exp(temp_a / temp_b);
end
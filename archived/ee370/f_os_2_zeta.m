function zeta_out = f_os_2_zeta(input_os)
    temp1 = -1 * log(input_os);
    temp2 = sqrt(pi^2 + (log(input_os))^2);
    zeta_out = temp1 / temp2;
end
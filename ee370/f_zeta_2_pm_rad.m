function pm_out = f_zeta_2_pm_rad(zeta_in)
    temp_a = 2 * zeta_in;
    temp_b = -2 * zeta_in^2;
    temp_c = sqrt(1 + 4 * zeta_in^4);
    temp_d = sqrt(temp_b + temp_c);
    pm_out = atan(temp_a / temp_d);
end
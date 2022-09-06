function to_m = f_cmil2m(c_mils)
    to_inches = c_mils / 1000;
    to_cm = to_inches * 2.54;
    to_m = 100 * to_cm;
end
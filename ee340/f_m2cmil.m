function to_cmil = f_m2cmil(meterz)
    to_cm = meterz / 100;
    to_inch = to_cm / 2.54;
    to_cmil = 1000 * to_inch;
end
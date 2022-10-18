function r_load = rload_from_Gamma(Gamma)
    tempn = 1 - real(Gamma)^2 - imag(Gamma)^2;
    tempd = (1 - real(Gamma))^2 + imag(Gamma)^2;
    r_load = tempn / tempd; 
end
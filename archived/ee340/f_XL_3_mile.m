%{
    3-phase @ 60 Hz
    reactance in line,   ohms / mile
%}
function reactance = f_XL_3_mile(GMD, GMR)
    reactance = 0.1213 * log(GMD/GMR);
end
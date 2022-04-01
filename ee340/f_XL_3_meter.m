%{
    3-phase @ 60 Hz
    reactance in line,   ohms / meter
%}
function reactance = f_XL_3_meter(GMD, GMR)
    reactance = (0.754 / 10^4) * log(GMD/GMR);
end
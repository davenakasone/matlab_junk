function [temp_y1, temp_y2, temp_y3] = f_delta2y(delta_A, delta_B, delta_C)
%{
    transform delta(Z) to Y(Z)
%}  
    temp_den = delta_A + delta_B + delta_C;
    temp_y1 = (delta_B * delta_C) / temp_den;
    temp_y2 = (delta_A * delta_C) / temp_den;
    temp_y3 = (delta_A * delta_B) / temp_den;
end
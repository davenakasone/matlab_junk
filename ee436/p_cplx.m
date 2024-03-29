function p_cplx(note, val)
    
if imag(val) < 0
    fprintf("%s:  %0.3f  %0.3f   |%0.3f|angle(%0.2f deg or %0.2f rad)\n",...
       note,...
       real(val), ...
       imag(val), ...
       abs(val), ...
       rad2deg(angle(val)), ...
       angle(val));
else
    fprintf("%s:  %0.3f + %0.3f   |%0.3f|angle(%0.2f deg or %0.2f rad)\n",...
       note,...
       real(val), ...
       imag(val), ...
       abs(val), ...
       rad2deg(angle(val)), ...
       angle(val));
end

end
%{
    for a solid core 3-phase, GMR = 0.7788 r

%}
function val_rp = f_line_rp(r_in)
    val_rp = r_in * exp(-1/4);
end
function r_t = f_routh(aa, bb, cc, dd, col)
    r_t = -1 * (aa*dd - bb*cc) / col;
end
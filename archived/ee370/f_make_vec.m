function vec_out = f_make_vec(pt_start, pt_stop)
    out_y = imag(pt_stop) - imag(pt_start);
    out_x = real(pt_stop) - real(pt_start);
    vec_out = out_x + 1j * out_y;
end
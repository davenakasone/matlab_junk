function tf_out = f_Fs_2_tf(Fs_in)
    [fs_n, fs_d] = numden(Fs_in);
    Fs_n = sym2poly(fs_n);
    Fs_d = sym2poly(fs_d);
    tf_out = tf(Fs_n, Fs_d);
end
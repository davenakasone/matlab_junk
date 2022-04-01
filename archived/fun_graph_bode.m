function fun_graph_bode(numerator, denominator,freqs)
%{
    numerator:    [ zn , zn-1, ..., z0]       coefs for zeros
    denominator:  [ pn , pn-1, ..., p0]       coefs for poles
    freq [start, stop]                        what range you want to see
    evals [ pt1, pt2, ..., ptn ]              specific evaluation freqencies
%}

H = tf(numerator, denominator);
hand = bode(H, { freqs(1), freqs(2) } );
grid on;
view(2);


end


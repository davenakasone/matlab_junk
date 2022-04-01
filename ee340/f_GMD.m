%{
    geometric mean distance
    (d1 d2 d3)^(1/3)
%}
function g_m_d = f_GMD(d_1, d_2, d_3)
    g_m_d = (d_1 * d_2 * d_3) .^ (1/3);
end
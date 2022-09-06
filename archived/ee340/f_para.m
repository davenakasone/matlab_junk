function imp_para = f_para(imp_a, imp_b)
%{
    provide 2 impedances
    get single parallel value
%}
    tempN = imp_a * imp_b;
    tempD = imp_a + imp_b;
    imp_para = tempN / tempD;
end
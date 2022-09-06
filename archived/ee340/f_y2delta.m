function [del_1, del_2, del_3] = f_y2delta(impy_1, impy_2, impy_3)
%{
    convert Y to delta
%}
    temp_num = impy_1*impy_2 + impy_2*impy_3 + impy_3*impy_1;
    del_1 = temp_num / impy_1;
    del_2 = temp_num / impy_2;
    del_3 = temp_num / impy_3;
end
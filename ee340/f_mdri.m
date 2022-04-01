function f_mdri(str_name, com_num, eng_mult)
%{
    fprintf a complex number
        name of number
        complex number to format
        stopping string/formats
        multiple
%}
    temp_mag = abs(com_num) * eng_mult;
    temp_deg = rad2deg(angle(com_num));
    temp_real = real(com_num) * eng_mult;
    temp_imag = imag(com_num) * eng_mult;
    fprintf("%s  , |%s|=  %0.4f  @  %0.4f deg  ,  real()=  %0.4f  ,  imag()=  %0.4f\n",...
        str_name, str_name, temp_mag, temp_deg, temp_real, temp_imag);
    %{
    fprintf("%s  , |%s|=  %0.3f  @  %0.2f deg  ,  real()=  %0.2f  ,  imag()=  %0.2f    multiple> %d\n",...
        str_name, str_name, temp_mag, temp_deg, temp_real, temp_imag, eng_mult);
    %}
end
%{
    final review

    1  :  Routh, special case 1, 0 in a column, use epsilon...try preparation also
    2  :  Routh, special case 2, row becomes 0
    

%}
format compact;
clear;
close all;
clc;
select = 2;


%------------------------------------------------------------------------------------------
if (select == 1)
    % given char: s^5 + 2*s^4 + 3*s^3 + 6*s^2 + 5*s + 3    CHANGE
    syms s;
    syms ep;
    characteristic = (s^5 + 2*s^4 + 3*s^3 + 6*s^2 +5*s + 3);
    characteristic_p = sym2poly(characteristic)
    quick_check = roots(characteristic_p)

    a_5 = 1; a_3 = 3; a_1 = 5;   %  0
    a_4 = 2; a_2 = 6; a_0 = 3;   %  0
    % c_11     c_12    % 0       %  0
    % c_21     c_22    % 0
    % c_31     c_32    % 0

    % first row, s^3
    c_11 = f_routh(a_5, a_3, a_4, a_2, a_4) % it == 0, so make it epsilon
    c_11 = ep;
    c_12 = f_routh(a_5, a_1, a_4, a_0, a_4);

    % second row, s^2
    c_21 = f_routh(a_4, a_2, c_11, c_12, c_11);
    c_22 = f_routh(a_4, a_0, c_11, 0, c_11);

    % third row, s^1
    c_31 = f_routh(c_11, c_12, c_21, c_22, c_21);
    c_32 = f_routh(c_11, 0, c_21, 0, c_21);

    % fourth row
    c_41 = f_routh(c_21, c_22, c_31, c_32, c_31);

    epz = linspace(-5, 5, 100);
    figure;
    grid on;
    hold on;
    plot(epz, subs(c_21, ep, epz), "r-", LineWidth=2);
    plot(epz, subs(c_31, ep, epz), "b-", LineWidth=2);
    hold off;

end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    Cs = s^5 + 7*s^4 + 6*s^3 + 42*s^2 + 8*s + 56;
    check_roots = roots(sym2poly(Cs))
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
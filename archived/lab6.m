%{
    1  :  exp2
    2  :  exp3
    3  :  exp4
    4  :  exp4
%}
format compact;
close all;
clear all;
clc;
select = 4;


%-------------------------------------------------------------------------------------
if select == 1
    I_f = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    E_A = [27.72, 53.14, 74.92, 93.26, 106.43, 116.16, 123.98, 129.03, 134.05];
    siz_arr1 = size(I_f);
    siz_arr2 = size(E_A);
    if siz_arr1 ~= siz_arr2
        fprintf("\nbad>> siz1:  %0.3f  ,  siz2:  %0.3f\n", siz_arr1(2), siz_arr2(2));
    end
    my_siz = siz_arr2(2);
    fprintf("I_f (A)  |  E_A  (V)\n");
    fprintf("---------------------\n");
    for ii = 1:1:my_siz
        fprintf("  %0.2f   |  %0.2f\n", I_f(ii), E_A(ii));
    end
    
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("E_A vs I_f", 'fontsize', 18);
    xlabel("I_f (A)", 'fontsize', 14);
    ylabel("E_A (V)", 'fontsize', 14);
    plot(I_f, E_A, 'b-', 'linewidth', 3);
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 2
    I_f = [0.05 , 0.1  , 0.15 , 0.2  , 0.25];
    I_A = [0.135, 0.249, 0.363, 0.481, 0.612];
    siz_arr1 = size(I_f);
    siz_arr2 = size(I_A);
    if siz_arr1 ~= siz_arr2
        fprintf("\nbad>> siz1:  %0.3f  ,  siz2:  %0.3f\n", siz_arr1(2), siz_arr2(2));
    end
    my_siz = siz_arr2(2);
    fprintf("I_f (A)  |  I_A  (A)\n");
    fprintf("---------------------\n");
    for ii = 1:1:my_siz
        fprintf("  %0.2f   |  %0.2f\n", I_f(ii), I_A(ii));
    end
    
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("I_A vs I_f", 'fontsize', 18);
    xlabel("I_f (A)", 'fontsize', 14);
    ylabel("I_A (A)", 'fontsize', 14);
    plot(I_f, I_A, 'b-', 'linewidth', 3);
    hold off;
end


%-------------------------------------------------------------------------------------
if select == 3
    I_f = [0.67, 0.67  , 0.67  , 0.67  , 0.67  , 0.67  , 0.67  ];
    I_A = [0   , 0.95  , 0.19  , 0.28  , 0.38  , 0.55  , 0.62  ];
    E_A = [208 , 205.43, 202.38, 199.36, 196.54, 190.23, 187.76];
    R =   [2000, 1200  , 600   , 400   , 300   , 300   , 171   ];
    fprintf("I_f (A)  |  R (ohm)  |  I_A (A)  |  E_A (V)\n");
    fprintf("---------------------------------------------\n");
    for ii = 1:1:7
        fprintf("  %0.2f   |  %7.2f  |  %0.2f     | %0.2f\n",...
            I_f(ii), R(ii), I_A(ii), E_A(ii));
    end
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("I_A, I_f by load resistance", 'fontsize', 18);
    xlabel("R (ohm)", 'fontsize', 14);
    ylabel("(A)", 'fontsize', 14);
    plot(flip(R), flip(I_f), 'b-', 'linewidth', 3);
    plot(flip(R), flip(I_A), 'r-', 'linewidth', 3);
    %plot(R, E_A, 'g-', 'linewidth', 3);
    legend('I_f', 'I_A', 'best', 'fontsize', 12);
    hold off;
    %
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("E_A by load resistance", 'fontsize', 18);
    xlabel("R (ohm)", 'fontsize', 14);
    ylabel("(V)", 'fontsize', 14);
    plot(flip(R), flip(E_A), 'g-', 'linewidth', 3);
    %plot(R, E_A, 'g-', 'linewidth', 3);
    legend('E_A', 'best', 'fontsize', 12);
    hold off;
    %}
end


%-------------------------------------------------------------------------------------
if select == 4
    I_f = [0.67, 0.68  , 0.71  , 0.77  , 0.83];
    I_A = [0   , 0.11  , 0.22  , 0.42  , .61 ];
    E_A = [208 , 208   , 208   , 208   , 208 ];
    R =   [2000, 1200  , 600   , 300   , 200 ];
    fprintf("I_f (A)  |  R (ohm)  |  I_A (A)  |  E_A (V)\n");
    fprintf("---------------------------------------------\n");
    for ii = 1:1:4
        fprintf("  %0.2f   |  %7.2f  |  %0.2f     | %0.2f\n",...
            I_f(ii), R(ii), I_A(ii), E_A(ii));
    end
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("I_A, I_f by load resistance", 'fontsize', 18);
    xlabel("R (ohm)", 'fontsize', 14);
    ylabel("(A)", 'fontsize', 14);
    plot(flip(R), flip(I_f), 'b-', 'linewidth', 3);
    plot(flip(R), flip(I_A), 'r-', 'linewidth', 3);
    %plot(R, E_A, 'g-', 'linewidth', 3);
    legend('I_f', 'I_A', 'best', 'fontsize', 12);
    hold off;
    %
    %
    figure('Position', [20, 20, 800, 800]);
    hold on;
    grid on;
    title("E_A by load resistance", 'fontsize', 18);
    xlabel("R (ohm)", 'fontsize', 14);
    ylabel("(V)", 'fontsize', 14);
    plot(flip(R), flip(E_A), 'g-', 'linewidth', 3);
    %plot(R, E_A, 'g-', 'linewidth', 3);
    legend('E_A', 'best', 'fontsize', 12);
    hold off;
    %}
end


%-------------------------------------------------------------------------------------
if select == 99
    
end
%%%%%%%%~~~~~~~~~END> 

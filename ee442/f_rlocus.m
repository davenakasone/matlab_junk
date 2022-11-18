function f_rlocus(input_tf)

    sz = size(input_tf);
    your_position = [20, 20, 800, 800];
    
    figure('Position', your_position);
    hold on;
    grid on;
    %axis equal;
    %axis padded;
    if sz(2) == 1
        rlocus(input_tf(1), 'g');
        legend('sys1');
    elseif sz(2) == 2
        rlocus(input_tf(1), 'g', input_tf(2), 'b');
        legend('sys1', 'sys2', fontsize=20);
    else
        rlocus(input_tf(1), 'g', input_tf(2), 'b', input_tf(3), 'r');
        legend('sys1', 'sys2', 'sys3');
    end
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
end
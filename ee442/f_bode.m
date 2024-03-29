function f_bode(input_tf)

    sz = size(input_tf);
    your_position = [30, 30, 800, 800];
    
    figure('Position', your_position);
    hold on;
    %axis equal;
    %axis padded;
    if sz(2) == 1
        bode(input_tf(1,1), 'g-');
        legend('sys1', fontsize=20);
    elseif sz(2) == 2
        bode(input_tf(1,1), 'g-', input_tf(1,2), 'b-');
        legend('sys1', 'sys2', fontsize=20);
    else
        bode(input_tf(1,1), 'g-', input_tf(1,2), 'b-', input_tf(1,3), 'r--');
        legend('sys1', 'sys2', 'sys3', fontsize=20);
    end
    set(findall(gcf,'type','line'),'linewidth',2);
    grid on;
    hold off;
end
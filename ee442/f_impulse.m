function f_impulse(input_tf)

    sz = size(input_tf);
    your_position = [50, 50, 800, 800];
    
    figure('Position', your_position);
    hold on;
    grid on;
    %axis equal;
    %axis padded;
    if sz(2) == 1
        impulse(input_tf(1,1), "g-");
        legend('sys1', fontsize=20);
    elseif sz(2) == 2
        impulse(input_tf(1,1), "g-", input_tf(1,2), "b-");
        legend('sys1', 'sys2', fontsize=20);
    else
        impulse(input_tf(1,1), "g-", input_tf(1,2), "b-", input_tf(1,3), 'r--');
        legend('sys1', 'sys2', 'sys3', fontsize=20);
    end
    set(findall(gcf,'type','line'),'linewidth',2);
    set(groot,'defaultLineMarkerSize',10);
    hold off;
end
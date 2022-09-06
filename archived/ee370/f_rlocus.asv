function f_rlocus(input_tf)

    your_position = [20, 20, 800, 800];
    
    figure('Position', your_position);
    hold on;
    grid on;
    axis equal;
    %axis padded;
    rlocus(input_tf);
    set(findall(gcf,'type','line'),'linewidth',3);
    set(groot,'defaultLineMarkerSize',15);
    sgrid;
    hold off;
end
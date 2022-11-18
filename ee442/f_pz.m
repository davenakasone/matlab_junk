function f_pz(input_tf)

    sz = size(input_tf);
    your_position = [40, 40, 800, 800];
    
    figure('Position', your_position);
    hold on;
    grid on;
    %axis equal;
    %axis padded;
    %pzplot(input_tf);
    popt = pzoptions;
    %popt.Grid = 'on';
    
    if sz(2) == 1
        pzmap(input_tf(1,1), 'g', popt);
        legend('sys1', fontsize=20);
        fprintf("\n");
        your_TF1 = input_tf(1)
        [tf_n, tf_d] = tfdata(input_tf(1), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
    elseif sz(2) == 2
        pzmap(input_tf(1,1), 'g', input_tf(1,2), 'b', popt);
        legend('sys1', 'sys2', fontsize=20);
        fprintf("\n");
        your_TF1 = input_tf(1)
        [tf_n, tf_d] = tfdata(input_tf(1), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
        fprintf("\n");
        your_TF2 = input_tf(2)
        [tf_n, tf_d] = tfdata(input_tf(2), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
    else
        pzmap(input_tf(1,1), 'g', input_tf(1,2), 'b', input_tf(1,3), 'r', popt);
        legend('sys1', 'sys2', 'sys3', fontsize=20);
        fprintf("\n");
        your_TF1 = input_tf(1)
        [tf_n, tf_d] = tfdata(input_tf(1), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
        fprintf("\n");
        your_TF2 = input_tf(2)
        [tf_n, tf_d] = tfdata(input_tf(2), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
        fprintf("\n");
        your_TF3 = input_tf(3)
        [tf_n, tf_d] = tfdata(input_tf(3), 'v');
        [z, p, k] = tf2zp(tf_n, tf_d)
    end
   
    set(findall(gcf,'type','line'),'linewidth',2);
    set(groot,'defaultLineMarkerSize',8);
    sgrid;
    hold off;
end
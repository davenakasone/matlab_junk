function fun_conTabg(yourFun, alpha, beta, gama, inpR, pnts)
%{
    yourFun as a function of GLOBAL pt     gamma * X( alpha * t + beta)

    alpha = 1       nothing
    alpha = 0       abortion, don't do
    alpha < 0       time reversal    flip on y-axis
    abs(alpha) < 1  horizontal stretch
    abs(alpha) > 1  horizontal compress

    beta = 0        nothing
    beta < 0        delay
    beta > 0        advance
    
    gama = 1        nothing
    gama < 0        flip on x axis
    abs(gama) > 1   vertical stretch
    abs(gama) < 1   vertical compress
        
    inpR [ min, max]    range for inputs , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]

    global pt won't be altered...it just gets transfered to a temporary t
%}
outR = [0,0];
buf = 1;
dots = 100;
szChk = size(pnts);
global pt;

    t = linspace(inpR(1), inpR(2), dots);
    original = subs(yourFun, pt, t);
    
    altered = zeros(1,dots);   % gamma * X( alpha * t + beta)
    for k = 1:dots
        altered(k) = gama * subs(yourFun, pt, t(k)*alpha + beta);
        if altered(k) > outR(2)
            outR(2) = altered(k);
        end
        if altered(k) < outR(1)
            outR(1) = altered(k);
        end
    end

    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    axis equal;
    view(2); % 2 for 2D
    tiStr = sprintf(' %.2f  * X ( %.2f t  +  %.2f )', gama, alpha, beta);
    title(tiStr, 'fontsize', 16); 
    xlabel('t');
    ylabel('X(t)');   
    xlim([inpR(1)-buf, inpR(2)+buf]);
    ylim([outR(1)-buf, outR(2)+buf]);
    indpAx = linspace(inpR(1)-buf , inpR(2)+buf , dots);
    dpntAx = linspace(outR(1)-buf , outR(2)+buf , dots);
    plot(indpAx  , 0*indpAx, 'k', 'linewidth', 1);
    plot(0*dpntAx, dpntAx  , 'k', 'linewidth', 1);        
    plot(inpR(2)+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , outR(2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    
    plot(t, original, 'r' , 'LineWidth', 3);
    %scatter(t, altered);
    %line(t, altered , 'b:', 'LineWidth', 5);
    plot(t, altered , 'b--', 'LineWidth', 2);
    %{
    for k = 1:dots
        plot(t(k), altered(k), 'b.', 'markersize', 10);
    end
    %}
    
    if szChk(2) > 1 % at least one point to plot
        for k = 1:szChk(1)
            plot(pnts(k,1), pnts(k,2), 'g.', 'markersize', 20);
            if pnts(k,1) < inpR(1) || pnts(k,1) > inpR(2)...
                    || pnts(k,2) < outR(1) || pnts(k,2) > outR(2)
                fprintf('\n\tpoint %d  ( %.1f , %.1f ) is out of range\n',...
                    k, pnts(k,1), pnts(k,2));
            end
        end
    end
    
end




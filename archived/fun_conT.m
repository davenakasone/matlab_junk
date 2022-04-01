function fun_conT(yourFun, inpR, pnts)
%{
    yourFun as a function of GLOBAL pt
    inpR [ min, max]    range for inputs  , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]

    global pt won't be altered...it just gets transfered to a temporary t
%}
outR = [0,0];
buf = 1;
szChk = size(pnts);
dots = 10 * (inpR(2)-inpR(1)+1);
global pt;

    t = linspace(inpR(1), inpR(2), dots);
    tempFun = subs(yourFun, pt, t);
    for k = 1:dots
        if tempFun(k) > outR(2)
            outR(2) = tempFun(k);
        end
        if tempFun(k) < outR(1)
            outR(1) = tempFun(k);
        end
    end

    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    %axis equal;
    view(2); % 2 for 2D
    tiStr = sprintf('%s', yourFun);
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
    
    plot(t, tempFun, 'r','LineWidth', 3);
    
    if szChk(2) > 1 % at least one point to plot
        for k = 1:szChk(1)
            plot(pnts(k,1), pnts(k,2), 'b.', 'markersize', 20);
            if pnts(k,1) < inpR(1) || pnts(k,1) > inpR(2)...
                    || pnts(k,2) < outR(1) || pnts(k,2) > outR(2)
                fprintf('\n\tpoint %d  ( %.1f , %.1f ) is out of range\n',...
                    k, pnts(k,1), pnts(k,2));
            end
        end
    end
    
end


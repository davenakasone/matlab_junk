function fun_conT2(funA, inpA, funB, inpB, pnts)
%{
    funA, funB as a functions of GLOBAL pt
    inpA, inpB [ min, max]    range for inputs  , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]

    global pt won't be altered...it just gets transfered to a temporary t
%}
outR = [0,0];
inpR = [0,0];
buf = 1;
szChk = size(pnts);
global pt;

    if inpA(1) >= inpA(2) || inpB(1) >= inpB(2)
        fprintf('fuck you bro\n');
    end
    if inpA(1) < inpB(1)
        inpR(1) = inpA(1);
    else
        inpR(1) = inpB(1);
    end
    if inpA(2) > inpB(2)
        inpR(2) = inpA(2);
    else
        inpR(2) = inpB(2);
    end
    dots = 10 * (inpR(2)-inpR(1));

    tA = linspace(inpA(1), inpA(2), dots);
    tB = linspace(inpB(1), inpB(2), dots);
    tempA = subs(funA, pt, tA);
    tempB = subs(funB, pt, tB);
    for k = 1:dots
        if tempA(k) > outR(2)
            outR(2) = tempA(k);
        end
        if tempB(k) > outR(2)
            outR(2) = tempB(k);
        end
        if tempA(k) < outR(1)
            outR(1) = tempA(k);
        end
        if tempB(k) < outR(1)
            outR(1) = tempB(k);
        end
    end

    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    %axis equal;
    view(2); % 2 for 2D
    tiStr = sprintf('first is red, second is blue');
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
    
    plot(tA, tempA, 'r','LineWidth', 2);
    plot(tB, tempB, 'b--','LineWidth', 2);
    
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



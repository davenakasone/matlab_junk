function fun_disN2(funA, inpA, funB, inpB, pnts)
%{
    funA, funB as a function of GLOBAL in
    inpA [ min, max]  и  inpB [ min, max]   range for inputs , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]
    globals won't be mutated
%}
outR = [0,0];
inpR = [0,0];
buf = 1;
szChk = size(pnts);
global in;

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

    nA = inpA(1):1:inpA(2);
    nSzA = size(nA);
    nB = inpB(1):1:inpB(2);
    nSzB = size(nB);
    
    if nSzA(2) > nSzB(2)
        dots = 10*nSzA(2);
    else
        dots = 10*nSzB(2);
    end
    
    tempA = subs(funA, in, nA);
    tempB = subs(funB, in, nB);
    %fillerA = zeros(1, nSzA(2));
    %fillerB = zeros(1, nSzB(2));
    %hldA = linspace(inpA(1), inpA(2), dots);
    %hldB = linspace(inpB(1), inpB(2), dots);
    %riderA = subs(funA, i_n, hldA);
    %riderB = subs(funB, i_n, hldB);
    
    for k = 1:nSzA(2)
        if tempA(k) > outR(2)
            outR(2) = tempA(k);
        end
        if tempA(k) < outR(1)
            outR(1) = tempA(k);
        end
    end
    for k = 1:nSzB(2)
        if tempB(k) > outR(2)
            outR(2) = tempB(k);
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
    tiStr = sprintf('two discrete functions');
    title(tiStr, 'fontsize', 16); 
    xlabel('n');
    ylabel('X[n]');   
    xlim([inpR(1)-buf, inpR(2)+buf]);
    ylim([outR(1)-buf, outR(2)+buf]);
    indpAx = linspace(inpR(1)-buf , inpR(2)+buf , dots);
    dpntAx = linspace(outR(1)-buf , outR(2)+buf , dots);
    plot(indpAx  , 0*indpAx, 'k', 'linewidth', 1);
    plot(0*dpntAx, dpntAx  , 'k', 'linewidth', 1);        
    plot(inpR(2)+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , outR(2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    
    %plot(hldA, riderA, 'r--','LineWidth', 1);
    %plot(hldB, riderB, 'm--','LineWidth', 1);
    
    stmA = stem(nA, tempA);
    stmA.LineStyle = ':';
    stmA.LineWidth = 1;
    stmA.Marker = 'd'
    stmA.MarkerSize = 10;
    stmA.MarkerFaceColor = 'none';
    stmA.MarkerEdgeColor = 'b';
    stmA.Color = 'b';
    temA = stmA.BaseLine;
    temA.Visible = 'off';
    
    stmB = stem(nB, tempB);
    %stmB.LineStyle = '-';
    stmB.LineWidth = 2;
    stmB.Marker = 'x'
    stmB.MarkerSize = 10;
    stmB.MarkerFaceColor = 'r';
    stmB.MarkerEdgeColor = 'r';
    stmB.Color = 'r';
    stmB.LineStyle = ':';
    temB = stmB.BaseLine;
    temB.Visible = 'off';
    
    %{
    quivA = quiver(nA,fillerA,fillerA,tempA);
    quivA.AutoScale = 'off';
    quivA.ShowArrowHead = 'off';
    quivA.Marker = '.';
    quivA.MarkerSize = 15;
    quivA.MaxHeadSize = 20;
    %quivA.AutoScaleFactor = 1;
    quivA.Color = 'b';
    quivA.LineWidth = 1;
    
    quivB = quiver(nB,fillerB,fillerB,tempB);
    quivB.AutoScale = 'off';
    quivB.ShowArrowHead = 'off';
    quivB.Marker = '.';
    quivB.MarkerSize = 15;
    quivB.MaxHeadSize = 20;
    %quivB.AutoScaleFactor = 1;
    quivB.Color = 'g';
    quivB.LineWidth = 1;
    %}
    
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




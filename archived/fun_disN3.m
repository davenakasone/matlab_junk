function fun_disN3(funA, inpA, funB, inpB, funC, inpC, pnts)
%{
    funA, B, C  as function of GLOBAL in
    inpA,B,C [ min, max]   range for inputs , output range adjusted automatically
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
    tempIn = [inpA; inpB; inpC];
    inpR(1) = min(tempIn, [], 'all');
    inpR(2) = max(tempIn, [], 'all');
   
    nA = inpA(1):1:inpA(2);
    nSzA = size(nA);
    nB = inpB(1):1:inpB(2);
    nSzB = size(nB);
    nC = inpC(1):1:inpC(2);
    nSzC = size(nC);
    
    dots = 10 * (inpR(2)-inpR(1));
    
    tempA = subs(funA, in, nA);
    tempB = subs(funB, in, nB);
    tempC = subs(funC, in, nC);
    %fillerA = zeros(1, nSzA(2));
    %fillerB = zeros(1, nSzB(2));
    %fillerC = zeros(1, nSzC(2));
    %hldA = linspace(inpA(1), inpA(2), dots);
    %hldB = linspace(inpB(1), inpB(2), dots);
    %hldC = linspace(inpC(1), inpC(2), dots);
    %riderA = subs(funA, in, hldA);
    %riderB = subs(funB, in, hldB);
    %riderC = subs(funC, in, hldC);
    
    tempOut = [tempA; tempB; tempC];
    outR(1) = min(tempOut, [], 'all');
    outR(2) = max(tempOut, [], 'all');

    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    %axis equal;
    view(2); % 2 for 2D
    tiStr = sprintf('three discrete functions');
    title(tiStr, 'fontsize', 16); 
    xlabel('n');
    ylabel('x[n]');   
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
    %plot(hldC, riderC, 'c--','LineWidth', 1);
    
    stmA = stem(nA, tempA);
    stmA.LineStyle = ':';
    set(stmA, 'Marker', 'o', 'MarkerSize', 15);
    stmA.LineWidth = 2;
    stmA.MarkerFaceColor = 'none';
    stmA.MarkerEdgeColor = 'b';
    stmA.Color = 'b';
    temA = stmA.BaseLine;
    temA.Visible = 'off';
    
    stmB = stem(nB, tempB);
    stmB.LineStyle = '-.';
    set(stmB, 'Marker', 's', 'MarkerSize', 10);
    stmB.LineWidth = 2;
    stmB.MarkerFaceColor = 'none';
    stmB.MarkerEdgeColor = 'r';
    stmB.Color = 'r';
    stmB.LineStyle = ':';
    temB = stmB.BaseLine;
    temB.Visible = 'off';
    
    stmC = stem(nC, tempC);
    stmC.LineStyle = '-';
    set(stmC, 'Marker', '^', 'MarkerSize', 10);
    stmC.LineWidth = 2;
    stmC.MarkerFaceColor = 'none';
    stmC.MarkerEdgeColor = 'g';
    stmC.Color = 'g';
    stmC.LineStyle = ':';
    temC = stmC.BaseLine;
    temC.Visible = 'off';
    
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
            plot(pnts(k,1), pnts(k,2), 'k.', 'markersize', 20);
            if pnts(k,1) < inpR(1) || pnts(k,1) > inpR(2)...
                    || pnts(k,2) < outR(1) || pnts(k,2) > outR(2)
                fprintf('\n\tpoint %d  ( %.1f , %.1f ) is out of range\n',...
                    k, pnts(k,1), pnts(k,2));
            end
        end
    end
    
end





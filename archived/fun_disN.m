function fun_disN(yourFun, inpR, pnts)
%{
    yourFun as a function of GLOBAL in
    inpR [ min, max]    range for inputs , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]

    global i_n won't be altered...it just gets transfered to a temporary n
%}
outR = [0,0];
buf = 1;
szChk = size(pnts);
global n;
global pt;

    nn = inpR(1):1:inpR(2);
    nSz = size(nn);
    dots = 10*nSz(2);
    tempFun = subs(yourFun, n, nn);
    filler = zeros(1, nSz(2));
    hld = linspace(inpR(1), inpR(2), dots);
    rider = subs(yourFun, n, hld);
    
    for k = 1:nSz(2)
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
    tiStr = sprintf('%s  n[ %.1f , %.1f ] w/ %.0f points',...
        yourFun, inpR(1), inpR(2), nSz(2) );
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
    
    %plot(hld, rider, 'g--','LineWidth', 1);
    
    %{
    quivH = quiver(n,filler,filler,tempFun);
    quivH.AutoScale = 'off';
    quivH.ShowArrowHead = 'off';
    quivH.Marker = '.';
    quivH.MarkerSize = 15;
    quivH.MaxHeadSize = 20;
    %quivH.AutoScaleFactor = 1;
    quivH.Color = 'b';
    quivH.LineWidth = 1;
    %}
    
    stmH = stem(nn, tempFun);
    set(stmH, 'Marker', 'o', 'MarkerSize', 10);
    stmH.LineStyle = '-';
    stmH.LineWidth = 2;
    stmH.MarkerFaceColor = 'none';
    stmH.MarkerEdgeColor = 'r';
    stmH.Color = 'b';
    tem = stmH.BaseLine;
    tem.Visible = 'off';
    
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




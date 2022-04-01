function fun_disNxxx(funs, range)
%{
    funs as function of GLOBAL in  [ fun1, fun2, fun3,....]
    range [ min, max]    for inputs , output range adjusted automatically
    globals won't be mutated, all functions get their own graph  use disNx if all on same graph
%}
buf = 1;
bump = 30;
pos = [20, 600, 700, 200];
points = 128;
szChk = size(funs);
outR = zeros(szChk(1,2),2);
global in;
    n = range(1):1:range(2);
    intv = range(2)-range(1)+1;
    tempFun = zeros(szChk(1,2),intv);
    for k = 1:szChk(1,2)
        for m = 1:intv
            tempFun(k,m) = subs(funs(1,k), in, n(1,m));
            if tempFun(k,m) < outR(k,1)
                outR(k,1) = tempFun(k,m);
            end
            if tempFun(k,m) > outR(k,2)
                outR(k,2) = tempFun(k,m);
            end
        end
    end
    
    indAx = linspace(range(1,1)-buf, range(1,2)+buf, points);
    for idx = szChk(1,2):-1:1
        tis = sprintf('%s', funs(1,idx));
        depAx = linspace(outR(idx,1)-buf, outR(idx,2)+buf, points);
        
        figure('Name', sprintf('function %d', idx),...
               'Position', pos,...
               'NumberTitle', 'off');
        hold on;
        grid on;
        view(2);
        xlabel('n or k');
        ylabel('x[n], y[n], h[n]');
        
        xlim([range(1)-buf, range(2)+buf]);
        ylim([outR(idx,1)-buf, outR(idx,2)+buf]);
        plot(indAx       , 0*indAx        , 'k' , 'linewidth' , 1);
        plot(0*depAx     , depAx          , 'k' , 'linewidth' , 1);        
        plot(range(2)+buf, 0              , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
        plot(0           , outR(idx,2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
        
        stmH = stem(n, tempFun(idx,:));
        stmH.LineStyle = '-';
        set(stmH, 'Marker', '^', 'MarkerSize', 5);
        stmH.LineWidth = 2;
        stmH.MarkerFaceColor = 'b';
        stmH.MarkerEdgeColor = 'b';
        stmH.Color = 'r';
        temH = stmH.BaseLine;
        temH.Visible = 'off';
        
        pos(1,1) = pos(1,1)+bump;
        pos(1,2) = pos(1,2)-bump;
        hold off;
    end
end

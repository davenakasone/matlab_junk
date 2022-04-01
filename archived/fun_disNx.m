function fun_disNx(funs, range)
%{
    funs as function of GLOBAL in  [ fun1, fun2, fun3,....]
    range [ min, max]    for inputs , output range adjusted automatically
    any points you want   [ t, x(t) ] ... [ n x 2]
    globals won't be mutated
%}
buf = 1;
figure('Position',[10, 120, 1450, 900]);
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

    for k = 1:szChk(1,2)
        %figure(k);  
        subplot(szChk(1,2),1,k);
        hold on;
        grid on;
        view(2); 
        xlabel('n');
        ylabel('x[n]');   
        xlim([range(1)-buf, range(2)+buf]);
        ylim([outR(k,1)-buf, outR(k,2)+buf]);
        indpAx = linspace(range(1)-buf , range(2)+buf , 128);
        dpntAx = linspace(outR(k,1)-buf , outR(k,2)+buf , 128);
        plot(indpAx  , 0*indpAx, 'k', 'linewidth', 1);
        plot(0*dpntAx, dpntAx  , 'k', 'linewidth', 1);        
        plot(range(2)+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
        plot(0      , outR(k,2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
        
        stmH = stem(n, tempFun(k,:));
        stmH.LineStyle = '-';
        set(stmH, 'Marker', '^', 'MarkerSize', 5);
        stmH.LineWidth = 2;
        stmH.MarkerFaceColor = 'b';
        stmH.MarkerEdgeColor = 'b';
        stmH.Color = 'r';
        temH = stmH.BaseLine;
        temH.Visible = 'off';
    end
end

   



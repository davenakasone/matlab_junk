function fun_conTx(funs, range)
%{
    funs = [fun1, fun2, ...] a functions of GLOBAL pt
    range = [ min, max]    range for inputs  , output range adjusted automatically
    global pt won't be altered...it just gets transfered to a temporary t
%}
buf = 1;
bump = 30;
pos = [20, 600, 700, 200];
points = 10;
szChk = size(funs);
outR = zeros(szChk(1,2),2);
global pt;
dots = points * (range(2)-range(1));
t = linspace(range(1), range(2), dots);
tempFuns = zeros(szChk(1,2), 2);
    for k = 1:szChk(1,2)
        for m = 1:dots
            tempFuns(k,m) = subs(funs(1,k), pt, t(1,m));
            if tempFuns(k,m) < outR(k,1)
                outR(k,1) = tempFuns(k,m);
            end
            if tempFuns(k,m) > outR(k,2)
                outR(k,2) = tempFuns(k,m);
            end
        end
    end
    
    indAx = linspace(range(1,1)-buf, range(1,2)+buf, dots);
    for idx = szChk(1,2):-1:1
        %tis = sprintf('%s', funs(1,idx));
        depAx = linspace(outR(idx,1)-buf, outR(idx,2)+buf, dots);
        
        figure('Name', sprintf('function %d', idx),...
               'Position', pos,...
               'NumberTitle', 'off');
        hold on;
        grid on;
        view(2);
        xlabel('t or tau');
        ylabel('x(t), y(t), h(t)');
        
        xlim([range(1)-buf, range(2)+buf]);
        ylim([outR(idx,1)-buf, outR(idx,2)+buf]);
        plot(indAx       , 0*indAx        , 'k' , 'linewidth' , 1);
        plot(0*depAx     , depAx          , 'k' , 'linewidth' , 1);        
        plot(range(2)+buf, 0              , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
        plot(0           , outR(idx,2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
        plot(t, tempFuns(idx,:), 'r','LineWidth', 2);
        
        pos(1,1) = pos(1,1)+bump;
        pos(1,2) = pos(1,2)-bump;
        hold off;
    end
end


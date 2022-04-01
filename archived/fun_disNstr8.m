function fun_disNstr8(inpt, dpnt, title_nm)
%{
    inpt = [ fun1 n values;
             fun2 n values;
               ...        ];

    dpnt = [ fun1 values;
             fun2 values;
               ...     ];

    get x[n] , h[n] , y[n]  ect across whatever range you feed in
        no subs required, just graphs
update this later to use a cell array
%}
    buf = 1;
    pos = [20, 600, 700, 200];
    points = 128;

        xRng = [ min(inpt, [], 'all') , max(inpt, [], 'all') ];
        yRng = [ min(dpnt, [], 'all') , max(dpnt, [], 'all') ];
        indAx = linspace(xRng(1,1)-buf, xRng(1,2)+buf, points);
        depAx = linspace(yRng(1,1)-buf, yRng(1,2)+buf, points);
        
        figure('Position', pos, 'NumberTitle', 'off');
        hold on;
        grid on;
        view(2);
        xlabel('independent');
        ylabel('dependent');
        title(title_nm);
        xlim([xRng(1)-buf, xRng(2)+buf]);
        ylim([yRng(1)-buf, yRng(2)+buf]);
        plot(indAx       , 0*indAx        , 'k' , 'linewidth' , 1);
        plot(0*depAx     , depAx          , 'k' , 'linewidth' , 1);        
        plot(xRng(2)+buf, 0          , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
        plot(0          , yRng(2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
        
        stmH = stem(inpt, dpnt);
        stmH.LineStyle = '-';
        set(stmH, 'Marker', '^', 'MarkerSize', 5);
        stmH.LineWidth = 2;
        stmH.MarkerFaceColor = 'b';
        stmH.MarkerEdgeColor = 'b';
        stmH.Color = 'r';
        temH = stmH.BaseLine;
        temH.Visible = 'off';
end

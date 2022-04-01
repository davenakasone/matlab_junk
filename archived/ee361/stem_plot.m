function stem_plot(x_vals, y_vals)

    outx = [min(x_vals, [], 'all'), max(x_vals, [], 'all')];
    outy = [min(y_vals, [], 'all'), max(y_vals, [], 'all')];
    bufx = 1;
    bufy = 0.5;
    dots = 100;

    figure('Position',[20, 20, 700, 700]);
    hold on;
    grid on;
    %axis equal;
    view(2); % 2 for 2D
    tiStr = "p(x)";
    title(tiStr, 'fontsize', 16); 
    xlabel('n');
    ylabel('p(x)');   
    xlim([outx(1)-bufx, outx(2)+bufx]);
    ylim([outy(1)-bufy, outy(2)+bufy]);
    x_ax = linspace(outx(1)-bufx , outx(2)+bufx , dots);
    y_ax = linspace(outy(1)-bufx , outy(2)+bufx , dots);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);        
    plot(outx(2)+bufx, 0           , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0           , outy(2)+bufy, 'y.', 'markersize', 20, 'linewidth', 10); 
    
    
    stmH = stem(x_vals, y_vals);
    set(stmH, 'Marker', 'o', 'MarkerSize', 10);
    stmH.LineStyle = '-';
    stmH.LineWidth = 2;
    stmH.MarkerFaceColor = 'none';
    stmH.MarkerEdgeColor = 'r';
    stmH.Color = 'b';
    tem = stmH.BaseLine;
    tem.Visible = 'off';
end




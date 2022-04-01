function [n,y] = fun_disNconvo(x,h, span)
%{
    provide the input function  x[n] defined in global in
    provide the reaction function h[n]
    the span is a lose range, make sure your intveral encloses x   и  h
    get the n list and y[n] values back
%}
    dust = 2;
    buf = .5;
    dist = span(2)-span(1)+1;
    n = span(1):1:span(2);
    y = zeros(1, dist);
    global in;
    global N;

    for i = 1:dist
        y(1,i) = symsum( subs(x, in, N) * subs(h, in, n(i)-N), N, span(1)*dust,span(2)*dust);
    end

    figure('Position',[10, 20, 700, 700]);
    hold on;
    grid on;
    title('convolution');
    view(2); 
    xlabel('n');
    ylabel('y[n]');   
    xlim([span(1)-buf, span(2)+buf]);
    ylim([min(y)-buf, max(y)+buf]);
    indpAx = linspace(span(1)-buf , span(2)+buf , 128);
    dpntAx = linspace(min(y)-buf , max(y)+buf , 128);
    plot(indpAx  , 0*indpAx, 'k', 'linewidth', 1);
    plot(0*dpntAx, dpntAx  , 'k', 'linewidth', 1);        
    plot(span(2)+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , max(y)+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs

    stmH = stem(n, y);
    stmH.LineStyle = '-';
    set(stmH, 'Marker', '^', 'MarkerSize', 5);
    stmH.LineWidth = 2;
    stmH.MarkerFaceColor = 'b';
    stmH.MarkerEdgeColor = 'b';
    stmH.Color = 'r';
    temH = stmH.BaseLine;
    temH.Visible = 'off';
end


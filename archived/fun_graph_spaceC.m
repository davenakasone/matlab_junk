function fun_graph_spaceC(params, range)
%{
    given params [ x(t) , y(t), z(t) ]
    range [t min, t max]
    space curve gets plotted            [ x, y, 0] for complex
%}
figName = 'fig1';
buf = 1;
global pt;
syms t;
 points = 50;
temp = subs(params, pt, t);
  
        tiStr = sprintf('x(t): %s   , y(t):  %s  ,  z(t): %s      t (%.1f, %.1f )',...
            temp(1), temp(2), temp(3), range(1), range(2));
        stepS = ( range(2) - range(1) ) / 100; % step for animation

        
        maxPx = 10;
        maxPy = 10;
        maxPz = 10;

            figure('Name', figName, 'Position',[20, 20, 700, 700]);
            hold on;
            grid on;
            axis equal;
            view(125,30); % 3  also good
            title(tiStr, 'fontsize', 16); 
            xlabel('x axis');
            ylabel('y axis'); 
            zlabel('z axis');  
            xlim([-maxPx-buf, maxPx+buf]);
            ylim([-maxPy-buf, maxPy+buf]);
            zlim([-maxPz-buf, maxPz+buf]);
            xax = linspace(-maxPx-buf   , maxPx+buf   , points);
            yax = linspace(-maxPy-buf   , maxPy+buf   , points);
            zax = linspace(-maxPz-buf   , maxPz+buf   , points);
            plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
            plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
            plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);        
            plot3(maxPx+buf, 0 , 0,'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
            plot3(0, maxPy+buf, 0,'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
            plot3(0, 0, maxPz+buf,'y.', 'markersize', 20, 'linewidth', 10); % z / mag, imag, real

            xt = temp(1);
            yt = temp(2);
            zt = temp(3);
            spcC = fplot3(xt, yt, zt);
            spcC.TRange = range;
            spcC.LineWidth = 2;
            spcC.Color = 'r';
            %{
            for i = range(1):stepS:range(2)   % animation???
                spcC.XFunction = xt;
                spcC.YFunction = yt;
                spcC.ZFunction = zt;
                drawnow
            end
            %}
            start = subs(temp, t, range(1));
            stop = subs(temp, t, range(2));
            startS = sprintf('start ( %.2f, %.2f, %.2f )', start(1), start(2), start(3));
            stopS = sprintf('stop (%.2f, %.2f, %.2f )', stop(1), stop(2), stop(3));
            plot3(start(1), start(2), start(3), 'g.', 'markersize', 20);
            plot3(stop(1), stop(2), stop(3), 'g.', 'markersize', 20);
            text(start(1)+.2, start(2)+.2, start(3)+.2, startS, 'FontSize', 14);
            text(stop(1)+.2, stop(2)+.2, stop(3)+.2, stopS, 'FontSize', 14);
end


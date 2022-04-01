function  fun_graph2Dcomplex(pointsIn)
%{
    just puts some points on a plane   
        z = x + jy format                ...put the principal first
%}
buf = 5;
pts = 128;
szChk = size(pointsIn);

if szChk(1,1) < 1
    fprintf('bad input\n');
else
    maxX = 0;
    maxY = 0;
    for k = 1:szChk(1)
        if abs(real(pointsIn(k))) > maxX
            maxX = abs(real(pointsIn(k)));
        end
        if abs(imag(pointsIn(k))) > maxY
            maxY = abs(imag(pointsIn(k)));
        end
    end
    
    figure();
    hold on;
    grid on;
    view(2); % 2 for 2D
    title('your points', 'fontsize', 16); 
    xlabel('real');
    ylabel('imaginary');   
    xlim([-maxX-buf, maxX+buf]);
    ylim([-maxY-buf, maxY+buf]);
    rax = linspace(-maxX-buf   , maxX+buf   , pts);
    iax = linspace(-maxY-buf   , maxY+buf   , pts);
    plot(rax  , 0*rax, 'k', 'linewidth', 1);
    plot(0*iax, iax  ,'k', 'linewidth', 1);        
    plot(maxX+buf, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / real direction
    plot(0      , maxY+buf, 'y.', 'markersize', 20, 'linewidth', 10); % +y / imag direction
    
    plot(real(pointsIn(1)), imag(pointsIn(1)), 'b.', 'markersize', 35, 'linewidth', 10);
    if abs(imag(pointsIn(1))) < 0
        tiStr = sprintf(' principal: ( %.3f - j %.3f)', real(pointsIn(1)), abs(imag(pointsIn(1))));
    else
        tiStr = sprintf(' principal: ( %.3f + j %.3f)', real(pointsIn(1)), imag(pointsIn(1)));
    end
    text(real(pointsIn(1)), imag(pointsIn(1)), tiStr, 'FontSize', 16);
    
    for k = 2:szChk(1)
        plot(real(pointsIn(k)), imag(pointsIn(k)), 'g.', 'markersize', 15, 'linewidth', 10);
    end
end


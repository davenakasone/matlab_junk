function  fun_CT_zer_pol(zers,pols, roc, side)
    buf = 1;
    pts_roc = 6;
    stuff = 128;
    pos = [20,20,600,600];
    tiStr = 'signal';
    rngR = [0,0]; % sigma or real(s) horizontal axis
    rngC = [0,0]; % jw or imag(s) vertical axis)

    temp = size(zers);
    szZ = temp(1,2);
    temp = size(pols);
    szP = temp(1,2);

    for ct = 1:szZ
        if real(zers(1,ct)) < rngR(1)
            rngR(1) = real(zers(1,ct));
        end
        if real(zers(1,ct)) > rngR(2)
            rngR(2) = real(zers(1,ct));
        end
        if imag(zers(1,ct)) < rngC(1)
            rngC(1) = imag(zers(1,ct));
        end
        if imag(zers(1,ct)) > rngC(2)
            rngC(2) = imag(zers(1,ct));
        end
    end
    
    for ct = 1:szP
        if real(pols(1,ct)) < rngR(1)
            rngR(1) = real(pols(1,ct));
        end
        if real(pols(1,ct)) > rngR(2)
            rngR(2) = real(pols(1,ct));
        end
        if imag(pols(1,ct)) < rngC(1)
            rngC(1) = imag(pols(1,ct));
        end
        if imag(pols(1,ct)) > rngC(2)
            rngC(2) = imag(pols(1,ct));
        end
    end
    
    figure('Position', pos);
    title(tiStr, 'fontsize', 16);
    hold on;
    grid on;
    view(2);
    xlabel('real(s) = \sigma', 'FontSize', 14);
    ylabel('imag(s) = j\omega','FontSize', 14);   
    xlim([rngR(1)-buf, rngR(2)+buf]);
    ylim([rngC(1)-buf, rngC(2)+buf]);
    rax = linspace(rngR(1)-buf, rngR(2)+buf, stuff);
    iax = linspace(rngC(1)-buf, rngC(2)+buf, stuff);
    plot(rax  , 0*rax, 'k', 'linewidth', 1);
    plot(0*iax, iax  , 'k', 'linewidth', 1);        
    plot(rngR(2)+buf, 0          , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0          , rngC(2)+buf, 'y.', 'markersize', 20, 'linewidth', 10); 
    
    if roc ~=4
        for ct = 1:szZ
        plot(real(zers(1,ct)), imag(zers(1,ct)),...
            'bo', 'MarkerSize', 10,...
            'LineWidth', 2);
        fprintf('zero # %d  of %d  at  (%.3f, j %.3f)\n',...
            ct, szZ, real(zers(1,ct)), imag(zers(1,ct)) );
        end
        fprintf('\n');
        for ct = 1:szP
            plot(real(pols(1,ct)), imag(pols(1,ct)),...
                'rx', 'MarkerSize', 10,...
                'LineWidth', 2);
            fprintf('pole # %d  of %d  at  (%.3f, j %.3f)\n',...
                ct, szP, real(pols(1,ct)), imag(pols(1,ct)) );
        end
    end
        
    if roc == 1 % right-sided
        region = linspace(rngC(1)-buf+.1,rngC(2)+buf-.1, pts_roc);
        for ct = 1:pts_roc
            plot([side;rngR(2)+buf],[region(ct);region(ct)],'g--', 'LineWidth',2);
        end
        for ct = 1:stuff
            plot(side,iax(1,ct), 'g.')
        end
    end
    
    if roc == 2 % left-sided
        region = linspace(rngC(1)-buf+.1,rngC(2)+buf-.1, pts_roc);
        for ct = 1:pts_roc
            plot([rngR(1)-buf;side],[region(ct);region(ct)],'g--', 'LineWidth',2);
        end
        for ct = 1:stuff
            plot(side,iax(1,ct), 'g.')
        end
    end
    
    if roc == 3 % 2-sided, user defined intverval
        region = linspace(rngC(1)-buf+.1,rngC(2)+buf-.1, pts_roc);
        for ct = 1:pts_roc
            plot([side(1);side(2)],[region(ct);region(ct)],'g--', 'LineWidth',2);
        end
        for ct = 1:stuff
            plot(side(1),iax(1,ct), 'g.')
            plot(side(2),iax(1,ct), 'g.')
        end
    end
    
    if roc == 4 % finite, convg for all
        region = linspace(rngC(1)-buf+.1,rngC(2)+buf-.1, pts_roc);
        for ct = 1:pts_roc
            quiver(0+.1,region(ct),.9, 0,'g-', 'LineWidth',2);
            quiver(0-.1,region(ct),-.9, 0,'g-', 'LineWidth',2);
        end
    end
end


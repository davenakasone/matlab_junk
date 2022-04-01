function  fun_DT_zer_pol(zers, pols, roc, side)
%{
    you should only be working with real magnitudes "r"

    sub before call, don't feed symbols...only doubles

%}
    buf = 1;
    setter = .5;
    stuff = 128;
    pos = [20,20,600,600];
    tiStr = 'signal on Z plane w/ unit circle';
    len_z = length(zers);
    len_p = length(pols);
    rng_z = [999,-999];
    rng_p = [999,-999];
    
    syms th;
    pts_th = 49;
    ll_th = 0;
    ul_th = 2*sym(pi);
    var_th = linspace(ll_th, ul_th, pts_th);
    len_qv = 12;
    qv = [ 1*sym(pi)/7  , 2*sym(pi)/7  , 3*sym(pi)/7,...
            4*sym(pi)/7  , 5*sym(pi)/7  , 6*sym(pi)/7,...
            8*sym(pi)/7  , 9*sym(pi)/7  , 10*sym(pi)/7,...
            11*sym(pi)/7 , 12*sym(pi)/7 , 13*sym(pi)/7 ];
    
    for ct = 1:len_z
        if abs(zers(1,ct)) < rng_z(1)
            rng_z(1) = abs(zers(1,ct));
        end
        if abs(zers(1,ct)) > rng_z(2)
            rng_z(2) = abs(zers(1,ct));
        end
    end
    for ct = 1:len_p
        if abs(pols(1,ct)) < rng_p(1)
            rng_p(1) = abs(pols(1,ct));
        end
        if abs(pols(1,ct)) > rng_p(2)
            rng_p(2) = abs(pols(1,ct));
        end
    end
    rng(1) = min( [rng_z(1), rng_p(1)] );
    rng(2) = max( [rng_z(2), rng_p(2)] );
    
    rng_ax = [0,0];
    if rng(2) < 1
        rng_ax(2)= 1 + buf;
    else
        rng_ax(2) = rng(2) + buf;
    end
    xy_ax = linspace(-rng_ax(2), rng_ax(2), stuff);
    
    
    figure('Position', pos);
    title(tiStr, 'fontsize', 16);
    hold on;
    grid on;
    axis equal;
    view(2); 
    xlabel('real', 'FontSize', 14);
    ylabel('imag','FontSize', 14);   
    xlim([-rng_ax(2), rng_ax(2)]);
    ylim([-rng_ax(2), rng_ax(2)]);
    plot(xy_ax  , 0*xy_ax, 'k', 'linewidth', 1);
    plot(0*xy_ax, xy_ax  , 'k', 'linewidth', 1);        
    plot(rng_ax(2), 0          , 'y.', 'markersize', 20, 'linewidth', 10); 
    plot(0          , rng_ax(2), 'y.', 'markersize', 20, 'linewidth', 10); 
    
    unitx = cos(th);
    ux = subs(unitx, th, var_th);
    unity = sin(th);
    uy = subs(unity, th, var_th);
    plot(ux,uy, 'c-', 'LineWidth', .5);
    %text(1, setter, '1', 'FontSize', 10);
    
    % put the zeros on the right (positive)
    if len_z > 0
        for ct = 1:len_z
            plot(real(zers(1,ct)), imag(zers(1,ct)), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
            %strg = sprintf('|zz| = %.2f', abs(zers(1,ct)));
            strg = sprintf('%.2f , %.2f', real(zers(1,ct)), imag(zers(1,ct)));
            %text(real(zers(1,ct)), imag(zers(1,ct)), strg, 'FontSize', 10);
        end
    end
    % put the poles on the left (negative)
    if len_p > 0
        for ct = 1:len_p
            plot(real(pols(1,ct)),imag(pols(1,ct)), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
            %strg = sprintf('|zp| = %.2f', abs(zers(1,ct)));
            strg = sprintf('%.2f , %.2f', real(pols(1,ct)), imag(pols(1,ct)));
            %text(real(pols(1,ct)), imag(pols(1,ct)), strg, 'FontSize', 10);
        end
    end
    
    % right-sided , ROC outside of outer-most pole
    if roc == 1 
        outx = rng_p(2)*cos(th);
        ox = subs(outx, th, var_th);
        outy = rng_p(2)*sin(th);
        oy = subs(outy, th, var_th);
        plot(ox,oy, 'LineWidth', 2, 'LineStyle', '--', 'Color', [.651,.651,.651] );
        for ct = 1:len_qv
            quiver(rng_p(2)*cos(qv(ct)), rng_p(2)*sin(qv(ct)),...
                   buf*cos(qv(ct))     , buf*sin(qv(ct)),...
                   'g-',...
                   'AutoScaleFactor', .8,...
                   'LineWidth', 2);
        end        
    end
    
    % left-sided , ROC inside of inner-most pole
    if roc == 2
        inx = rng_p(1)*cos(th);
        ix = subs(inx, th, var_th);
        iny = rng_p(1)*sin(th);
        iy = subs(iny, th, var_th);
        plot(ix,iy, 'LineWidth', 2, 'LineStyle', '--', 'Color', [.651,.651,.651] );
        for ct = 1:len_qv
            quiver(rng_p(1)*cos(qv(ct)), rng_p(1)*sin(qv(ct)),...
                   -buf*cos(qv(ct))    , -buf*sin(qv(ct)),...
                   'g-',...
                   'AutoScaleFactor', 1,...
                   'LineWidth', 1);
        end        
    end
    
    % 2-sided, user defined intverval
    if roc == 3 
        inx = abs(side(1))*cos(th);
        ix = subs(inx, th, var_th);
        iny = abs(side(1))*sin(th);
        iy = subs(iny, th, var_th);
        plot(ix,iy, 'LineWidth', 2, 'LineStyle', '--', 'Color', [.651,.651,.651] );
        outx = abs(side(2))*cos(th);
        ox = subs(outx, th, var_th);
        outy = abs(side(2))*sin(th);
        oy = subs(outy, th, var_th);
        plot(ox,oy, 'LineWidth', 2, 'LineStyle', '--', 'Color', [.651,.651,.651] );
        for ct = 1:len_qv
            plot([subs(inx, th, qv(ct)), subs(outx, th, qv(ct))],...
                 [subs(iny, th, qv(ct)), subs(outy, th, qv(ct))],...
                 'LineWidth', 1,...
                 'LineStyle', '-',...
                 'Color', 'g');
        end
    end
end




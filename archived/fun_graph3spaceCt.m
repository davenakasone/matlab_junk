function    fun_graph3spaceCt(tipIn,range, points, sd)
%{
    provide "tipIn" as a vector that traces the path of a space curve
        [ x(t) , y(t) , z(t) ]      use global parameter t  "pt"

    provide "range" as start and stop of parameter

    any points you want plotted

    sd = 1 "static" just a few tip traces
    sd = 2 "dynamic" animation plays
    sd = 3  use slider at will
    sd = something else, only a few points get plotted

    don't try and make the viewpoint switch as time advances if you want a movie
    ....figH = figure();     flipBook(k) = getframe(figH, [20,20, 1000, 1000])  and it is good
%}
buf = 1;         % boundary to view graph
mult = 5;        % data points graphed
stuffing = 64;  % making axes
drop = 0;        % for label position
spin = .75;      % how many spins in view
pas = .2;        % how many seconds to pause
post = [20, 20, 1000, 1000];

szPoints = size(points);
rngX = [0,0];
rngY = [0,0];
rngZ = [0,0];
global pt;

if range(1) >= range(2)
    fprintf('\nthere is a problem with the range\n');
end
if (range(2) - range(1)) < 1
    dots = mult;
else
    dots = round( (range(2) - range(1))*mult);
end

    t = linspace(range(1), range(2), dots);
    scX = zeros(dots,1);
    scY = zeros(dots,1);
    scZ = zeros(dots,1);
    tX  = zeros(dots,1);
    tY  = zeros(dots,1);
    tZ  = zeros(dots,1);
    for k = 1:dots
        scX(k,1) = subs(tipIn(1), pt, t(k)+eps); % x(t)
        scY(k,1) = subs(tipIn(2), pt, t(k)+eps); % y(t)
        scZ(k,1) = subs(tipIn(3), pt, t(k)+eps); % z(t)
        tX(k,1) = subs( diff(tipIn(1), pt, 1), pt, t(k)+eps); % x'(t)
        tY(k,1) = subs( diff(tipIn(2), pt, 1), pt, t(k)+eps); % y'(t)
        tZ(k,1) = subs( diff(tipIn(3), pt, 1), pt, t(k)+eps); % z'(t)
        if rngX(1) > scX(k,1)
            rngX(1) = scX(k,1);
        end
        if rngX(2) < scX(k,1)
            rngX(2) = scX(k,1);
        end
        if rngY(1) > scY(k,1)
            rngY(1) = scY(k,1);
        end
        if rngY(2) < scY(k,1)
            rngY(2) = scY(k,1);
        end
        
        if rngZ(1) > scZ(k,1)
            rngZ(1) = scZ(k,1);
        end
        if rngZ(2) < scZ(k,1)
            rngZ(2) = scZ(k,1);
        end
    end
    start = [ scX(1,1)   , scY(1,1)   , scZ(1,1) ];
    stop =  [ scX(dots,1), scY(dots,1), scZ(dots,1) ];
    
%%                                  STATIC   sd == 1      
    if sd == 1
        tiStr = sprintf(' [ x(t) =  %s , y(t) =  %s , z(t) = %s ] position vector, t =  %.2f  to %.2f',...
        subs(tipIn(1), pt, 't'), subs(tipIn(2), pt, 't'), subs(tipIn(3), pt, 't'),...
        range(1), range(2));
        figure('Position', post);
        hold on;
        grid on;
        %axis equal;  % only if you need symmetry
            %view(90,0);   % yz
            %view(180,0);   %xz
            %view(0,90);    %  xy
            view(50,20); % 3  also good
        title(tiStr, 'fontsize', 16); 
        xlabel('x axis');
        ylabel('y axis'); 
        zlabel('z axis');  
        xlim([rngX(1)-buf, rngX(2)+buf]);
        ylim([rngY(1)-buf, rngY(2)+buf]);
        zlim([rngZ(1)-buf, rngZ(2)+buf]);
        xax = linspace(rngX(1)-buf, rngX(2)+buf , stuffing);
        yax = linspace(rngY(1)-buf, rngY(2)+buf , stuffing);
        zax = linspace(rngZ(1)-buf, rngZ(2)+buf , stuffing);
        plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
        plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        plot3(rngX(2)+buf, 0          , 0           ,'y.', 'markersize', 20, 'linewidth', 10); 
        plot3(0          , rngY(2)+buf, 0           ,'y.', 'markersize', 20, 'linewidth', 10);
        plot3(0          , 0          , rngZ(2)+buf ,'y.', 'markersize', 20, 'linewidth', 10);
        text(rngX(2)+buf-drop, drop            , drop            , '+ X');
        text(drop            , rngY(2)+buf-drop, drop            , '+ Y');
        text(drop            , drop            , rngZ(2)+buf-drop, '+ Z');
        %{
        zFill = zeros(dots,1);
        quivH = quiver3(zFill, zFill, zFill, scX, scY, scZ);
        quivH.AutoScale = 'on';
        quivH.AutoScaleFactor = 1;
        quivH.Color = 'b';
        quivH.LineWidth = .5;
        quivH.LineStyle = '-.';
        %quivH.Marker = '.';
        quivH.ShowArrowHead = 'off';
        %}
        %{
        quivF = quiver3(scX, scY, scZ, tX, tY, tZ);
        quivF.AutoScale = 'on';
        quivF.AutoScaleFactor = 1;
        quivF.Color = 'g';
        quivF.LineWidth = 2;
        quivF.LineStyle = '-';
        %quivF.Marker = '.';
        quivF.ShowArrowHead = 'on';
        %}
        plot3(scX, scY, scZ, 'r-', 'LineWidth', 2);
        
        plot3(start(1,1), start(1,2), start(1,3), 'k.', 'markersize', 5, 'linewidth', 10);
        text( start(1,1), start(1,2), start(1,3), 'start', 'FontSize', 16);
        plot3(stop(1,1) , stop(1,2) , stop(1,3) , 'k.', 'markersize', 5, 'linewidth', 10);
        text( stop(1,1) , stop(1,2) , stop(1,3) , 'stop', 'FontSize', 16);
        if szPoints(2) > 1 % at least one point to plot
            for k = 1:szPoints(1)
                if points(k,1) < rngX(1)-buf || points(k,1) > rngX(2)+buf ||...
                   points(k,2) < rngY(1)-buf || points(k,2) > rngY(2)+buf ||...
                   points(k,3) < rngZ(1)-buf || points(k,3) > rngZ(2)+buf
                       fprintf('\n\tpoint %d  ( %.1f , %.1f , %.1f ) is out of range\n',...
                           k, points(k,1), points(k,2), points(k,3) );
                else
                    plot3(points(k,1), points(k,2), points(k,3), 'g.', 'markersize', 20);
                end
            end
        end
        %
        % for hw3
        A = 20;
        %densX = linspace(
        %temp_x = linspace(rngX(1), rngX(2), dots);
            %temp_y = linspace( rngY(1), rngY(2), dots);
            %temp_z = linspace( rngZ(1), rngZ(2), dots);
            %dens_x = zeros(dots, dots) + A;
            %[dens_y, dens_z] = meshgrid(temp_y, temp_z);
        %[dens_x, dens_y, dens_z] = meshgrid(temp_x, temp_y, temp_z);
        %dens_x = A .* dens_x;
        %surf(dens_x, dens_y, dens_z,...
        %    'FaceColor', 'g',...
        %    'FaceAlpha', .3,...
        %    'EdgeColor', 'none'); % value of density function
        %
    end
    
    
%%                          DYNAMIC  sd == 2    
    if sd == 2
        figH = figure('Position', post);
        for k = 1:dots
            clf;
            tiStr = sprintf(' t = %.2f', k);
            hold on;
            grid on;
            %axis equal;  % only if you need symmetry
                %view(3);
                view(45,15); % use this or 3 if you want it static...turn off spin
                %view(spin*k*360/dots,35);
                %view(k*360/dots,k*90/dots);
                %view(spin, 30);
                %spin = spin + spin;
            title(tiStr, 'fontsize', 16); 
            xlabel('x axis');
            ylabel('y axis'); 
            zlabel('z axis');  
            xlim([rngX(1)-buf, rngX(2)+buf]);
            ylim([rngY(1)-buf, rngY(2)+buf]);
            zlim([rngZ(1)-buf, rngZ(2)+buf]);
            xax = linspace(rngX(1)-buf, rngX(2)+buf , stuffing);
            yax = linspace(rngY(1)-buf, rngY(2)+buf , stuffing);
            zax = linspace(rngZ(1)-buf, rngZ(2)+buf , stuffing);
            plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
            plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
            plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
            plot3(rngX(2)+buf, 0          , 0           ,'y.', 'markersize', 20, 'linewidth', 10); 
            plot3(0          , rngY(2)+buf, 0           ,'y.', 'markersize', 20, 'linewidth', 10);
            plot3(0          , 0          , rngZ(2)+buf ,'y.', 'markersize', 20, 'linewidth', 10);
            text(rngX(2)+buf-drop, drop            , drop            , '+ X');
            text(drop            , rngY(2)+buf-drop, drop            , '+ Y');
            text(drop            , drop            , rngZ(2)+buf-drop, '+ Z');
            
            plot3(start(1,1), start(1,2), start(1,3), 'k.', 'markersize', 5, 'linewidth', 10);
            text( start(1,1), start(1,2), start(1,3), 'start', 'FontSize', 16);
            plot3(stop(1,1) , stop(1,2) , stop(1,3) , 'k.', 'markersize', 5, 'linewidth', 10);
            text( stop(1,1) , stop(1,2) , stop(1,3) , 'stop', 'FontSize', 16);
            if szPoints(2) > 1 % at least one point to plot
                for m = 1:szPoints(1)
                    if points(m,1) < rngX(1)-buf || points(m,1) > rngX(2)+buf ||...
                       points(m,2) < rngY(1)-buf || points(m,2) > rngY(2)+buf ||...
                       points(m,3) < rngZ(1)-buf || points(m,3) > rngZ(2)+buf
                           fprintf('\n\tpoint %d  ( %.1f , %.1f , %.1f ) is out of range\n',...
                               m, points(m,1), points(m,2), points(m,3) );
                    else
                        plot3(points(m,1), points(m,2), points(m,3), 'g.', 'markersize', 20);
                    end
                end
            end
            
            plot3(scX, scY, scZ, 'r-', 'LineWidth', 2); % trajectory
            plot3(scX(k,1), scY(k,1), scZ(k,1), 'c.', 'markersize', 40); % current point
            quivH = quiver3(0,0,0,scX(k,1), scY(k,1), scZ(k,1));
            quivH.AutoScale = 'on';
            quivH.AutoScaleFactor = 1;
            quivH.Color = 'b';
            quivH.LineWidth = 1;
            quivH.LineStyle = '-';
            %quivH.Marker = '.';
            quivH.ShowArrowHead = 'off';
            %
            quivF = quiver3(scX(k,1), scY(k,1), scZ(k,1), tX(k,1), tY(k,1), tZ(k,1));
            quivF.AutoScale = 'on';
            quivF.AutoScaleFactor = 1;
            quivF.Color = 'g';
            quivF.LineWidth = 2;
            quivF.LineStyle = '-';
            %quivF.Marker = '.';
            quivF.ShowArrowHead = 'on';
            %}
            plot3(scX, scY, scZ, 'r-', 'LineWidth', 2);
                %drawnow;    % use these if you don't want the movie   
                pause(pas);      % use both if matlab isn't giving drawnow good behavior 
                %flipBook(k) = getframe(figH);
        end
    %{   
    vid = VideoWriter('spaceC', 'MPEG-4');
    vid.FrameRate = 3;
    open(vid);
    writeVideo(vid, flipBook);
    close(vid);
    %}
    end
    
    %%          SLIDER  sd == 3
    if sd == 3
        count = 3;
        xax = linspace(rngX(1)-buf, rngX(2)+buf , stuffing);
        yax = linspace(rngY(1)-buf, rngY(2)+buf , stuffing);
        zax = linspace(rngZ(1)-buf, rngZ(2)+buf , stuffing);
        figH = figure('Position', post);
        view(55,15);
        aes = gca;
        set(aes, 'xlim', [rngX(1)-buf, rngX(2)+buf]);
        set(aes, 'ylim', [rngY(1)-buf, rngY(2)+buf]);
        set(aes, 'zlim', [rngZ(1)-buf, rngZ(2)+buf]);
        xlabel('x axis');
        ylabel('y axis'); 
        zlabel('z axis');
        view(aes, 3);
        hold(aes, 'on');
        grid(aes, 'on');
        plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
        plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
        plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
        plot3(rngX(2)+buf, 0          , 0           ,'y.', 'markersize', 20, 'linewidth', 10); 
        plot3(0          , rngY(2)+buf, 0           ,'y.', 'markersize', 20, 'linewidth', 10);
        plot3(0          , 0          , rngZ(2)+buf ,'y.', 'markersize', 20, 'linewidth', 10);
        text(rngX(2)+buf-drop, drop            , drop            , '+ X');
        text(drop            , rngY(2)+buf-drop, drop            , '+ Y');
        text(drop            , drop            , rngZ(2)+buf-drop, '+ Z');

        plot3(start(1,1), start(1,2), start(1,3), 'k.', 'markersize', 5, 'linewidth', 10);
        text( start(1,1), start(1,2), start(1,3), 'start', 'FontSize', 16);
        plot3(stop(1,1) , stop(1,2) , stop(1,3) , 'k.', 'markersize', 5, 'linewidth', 10);
        text( stop(1,1) , stop(1,2) , stop(1,3) , 'stop', 'FontSize', 16);
        if szPoints(2) > 1 % at least one point to plot
            for m = 1:szPoints(1)
                if points(m,1) < rngX(1)-buf || points(m,1) > rngX(2)+buf ||...
                   points(m,2) < rngY(1)-buf || points(m,2) > rngY(2)+buf ||...
                   points(m,3) < rngZ(1)-buf || points(m,3) > rngZ(2)+buf
                       fprintf('\n\tpoint %d  ( %.1f , %.1f , %.1f ) is out of range\n',...
                           m, points(m,1), points(m,2), points(m,3) );
                else
                    plot3(points(m,1), points(m,2), points(m,3), 'g.', 'markersize', 20);
                end
            end
        end
        traj = plot3(scX, scY, scZ, 'r-', 'LineWidth', 3);

        sld = uicontrol('Parent',figH,...
                        'Style','slider',...
                        'Position',[post(1)+40,post(2)+20,post(3)/3,10],...
                        'value',t(1),...
                        'min', t(1),...
                        'max', t(dots),...
                        'BackgroundColor', figH.Color );

        sldMin = uicontrol('Parent',figH,...
                           'Style','text',...
                           'Position',[post(1)+20,post(2)+5,30,30],...
                           'String',sprintf('%.1f',t(1)),...
                           'BackgroundColor',figH.Color);

        sldMax = uicontrol('Parent',figH,...
                           'Style','text',...
                           'Position',[40+post(1)+post(3)/3, post(2)+5,40,20],...
                           'String',sprintf('%.1f',t(dots)),...
                           'BackgroundColor',figH.Color);

        sldTil = uicontrol('Parent',figH,...
                           'Style','text',...
                           'Position',[(post(1)+(post(3)/3))/2,post(2)+1,15,15],...
                           'String','t',...
                           'BackgroundColor',figH.Color);

        plane_pts = 5;
        store = 1000;
            %px = zeros(plane_pts, plane_pts);
        py = zeros(plane_pts, plane_pts);
        x = linspace(-store, store, plane_pts);
            %y = linspace(-store, store, plane_pts);
        z = linspace(-store, store, plane_pts);
        [px, pz] = meshgrid(x, z);
            %[py, pz] = meshgrid(y, z);
        while count ~= 2
            k = sld.Value;              
            tiStr = sprintf(' t = %.2f', k);
            title(tiStr, 'fontsize', 16);
            
            pos = [subs(tipIn(1),pt,k), subs(tipIn(2),pt,k), subs(tipIn(3),pt,k)];
            cp = plot3(pos(1), pos(2), pos(3), 'c.', 'markersize', 40); % current point
            R = quiver3(0,0,0,pos(1), pos(2), pos(3)); % position vector
            R.AutoScale = 'on';
            R.AutoScaleFactor = 1;
            R.Color = 'b';
            R.LineWidth = 1;
            R.LineStyle = '-';
            %R.Marker = '.';
            R.ShowArrowHead = 'off';
            
            tan = [ subs(diff(tipIn(1), pt, 1), pt, k),...
                    subs(diff(tipIn(2), pt, 1), pt, k),...
                    subs(diff(tipIn(3), pt, 1), pt, k) ];
            dR = quiver3(pos(1), pos(2), pos(3),tan(1), tan(2), tan(3));
            dR.AutoScaleFactor = 1;
            dR.Color = 'g';
            dR.LineWidth = 2;
            dR.LineStyle = '-';
            %dR.Marker = '.';
            dR.ShowArrowHead = 'on';
            
            
            for rows = 1:plane_pts
                for cols = 1:plane_pts
                        %px(rows, cols) = 3*pos(1) + (2*pos(1)^3) - 2*pos(1)*py(rows,cols);
                    py(rows, cols) = 1.5 + pos(1)^2 - px(rows, cols)/(2*pos(1));
                end
            end
                    
            plnH = surf(px, py, pz, 'FaceColor', 'y', 'FaceAlpha', .3, 'EdgeColor', 'none'); 
        
            %{
            nor = cross(pos,tan)./norm(cross(pos,tan));
            nH = quiver3(pos(1), pos(2), pos(3),nor(1), nor(2), nor(3));
            nH.AutoScaleFactor = 1;
            nH.Color = 'm';
            nH.LineWidth = 2;
            nH.LineStyle = '-';
            %nH.Marker = '.';
            nH.ShowArrowHead = 'on';
            %}
            
            drawnow;
            %pause(pas);
            delete(cp); delete(R); delete(dR); delete(plnH); %delete(plnP); %delete(nH); 
        end
    end
end


function out = burn_fun(var)

    out = var + 1;


end











%{

%%  preliminary work
buf = 1;         % boundary to view graph
mult = 3;        % data points graphed
stuffing = 256;  % making axes
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
%% graphing
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
    traj = plot3(scX, scY, scZ, 'r-', 'LineWidth', 2);
    
    sld = uicontrol('Parent',figH,...
                    'Style','slider',...
                    'Position',[post(1)+20,post(2)+20,post(3)/3,10],...
                    'value',t(1),...
                    'min', t(1),...
                    'max', t(dots),...
                    'BackgroundColor', figH.Color );
  
    sldMin = uicontrol('Parent',figH,...
                       'Style','text',...
                       'Position',[post(1)+20,post(2)+10,20,20],...
                       'String',sprintf('%.1f',t(1)),...
                       'BackgroundColor',figH.Color);
                   
    sldMax = uicontrol('Parent',figH,...
                       'Style','text',...
                       'Position',[post(1)+post(3)/3, post(2)+10,20,20],...
                       'String',sprintf('%.1f',t(dots)),...
                       'BackgroundColor',figH.Color);
                   
    sldTil = uicontrol('Parent',figH,...
                       'Style','text',...
                       'Position',[(post(1)+(post(3)/3))/2,post(2)+1,15,15],...
                       'String','t',...
                       'BackgroundColor',figH.Color);
                
    while count ~= 2
    k = sld.Value;              
    tiStr = sprintf(' t = %.2f', k);
    title(tiStr, 'fontsize', 16);
    cp = plot3(subs(tipIn(1),pt,k), subs(tipIn(2),pt,k), subs(tipIn(3),pt,k),...
        'c.', 'markersize', 40); % current point
    R = quiver3(0,0,0,subs(tipIn(1),pt,k), subs(tipIn(2),pt,k), subs(tipIn(3),pt,k)); % position
    R.AutoScale = 'on';
    R.AutoScaleFactor = 1;
    R.Color = 'b';
    R.LineWidth = 1;
    R.LineStyle = '-';
    %R.Marker = '.';
    R.ShowArrowHead = 'off';
    
    dR = quiver3(subs(tipIn(1),pt,k), subs(tipIn(2),pt,k), subs(tipIn(3),pt,k),...
        subs(diff(tipIn(1), pt, 1), pt, k), subs(diff(tipIn(2), pt, 1), pt, k),...
        subs(diff(tipIn(3), pt, 1), pt, k)   );
    dR.AutoScale = 'on';
    dR.AutoScaleFactor = 1;
    dR.Color = 'g';
    dR.LineWidth = 2;
    dR.LineStyle = '-';
    %dR.Marker = '.';
    dR.ShowArrowHead = 'on';
        
    drawnow;
    pause(pas);
    delete(cp); delete(R); delete(dR);
    end

x=0:1e-2:2*pi;
y=sin(x);
dx=2;

a=gca;
p=plot(x,y);

set(gcf,'doublebuffer','on');

set(a,'xlim',[0 dx]);
set(a,'ylim',[min(y) max(y)]);

pos=get(a,'position');
Newpos=[pos(1), pos(2)-0.1, pos(3), 0.05];

xmax=max(x);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
% Setting up callback string to modify XLim of axis (gca)
% based on the position of the slider (gcbo)

h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-dx);

%}

%{
zeta = .5;                           % Damping Ratio
wn = 2;                              % Natural Frequency
sys = tf(wn^2,[1,2*zeta*wn,wn^2]); 
f = figure;

ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
h = stepplot(ax,sys);
setoptions(h,'XLim',[0,10],'YLim',[0,2]);

b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',zeta, 'min',0, 'max',1);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','Damping Ratio','BackgroundColor',bgcolor);
            
b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2])); 
%}
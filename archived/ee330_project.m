%{
    ee330 project (reusable fuse)

    #1: steady state, open
    #2: steady state, short  (fuse was activated, break-down)
    #3: transient --> Laplace transform
    #4: display the fuse/capacitor
    #5: angle, distance optimization
    #6: capacitance with some
    #7: component selection
    #8: final fuze

%}
clc;
close all;
clearvars;
sympref('MatrixWithSquareBrackets', 1);
sympref('PolynomialDisplayStyle', 'ascend');
old_val = sympref('HeavisideAtOrigin', 1);

%plot: 4, 5, 8
                sel = 4;  % CHANGE CHANGE CHANGE
        
global sig; syms sig; assume(sig, 'real'); % resistivity , sub in and out with const
global c1; syms c1; % just a holder for integration constant, solve for latter
global c2; syms c2; % just a holder for integration constant, solve for latter
global c3; syms c3; % just a holder for integration constant, solve for latter
global c4; syms c4; % just a holder for integration constant, solve for latter
global rx; global ry; global rz; % "rectangular x"   "rectangular y"  "rectangular z"       
global cr; global cf; global cz; % "cylindrical rho" "cylindrical fi" "cylindrical z"       fi = phi 
global sr; global st; global sf; % "spherical radius" "spherical theta" "spherical fi"   to avoid confusion
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); assume(cr >= 0);
syms cf; assume(cf, 'real'); assume(cf >= 0);
syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); assume(sr >= 0);
syms st; assume(st, 'real'); assume(st >= 0);
syms sf; assume(sf, 'real'); assume(sf >= 0);
global pL; syms pL; assume(pL, 'real'); assume(pL >= 0); % line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); assume(pS >= 0); % surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); assume(pV >= 0); % volume charge denisty in C / m^3
global ep0; syms ep0; assume(ep0, 'real'); % sub out with const.ep0 at end of calculation (if needed)
global chie; syms chie; assume(chie, 'real'); % as in  chie = epr - 1   or chie = ep/ep0 - 1
global pt; syms pt; assume(pt, 'real'); % paramater t for a space curve or "time" var
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect

%{
    ensure d > 0 so there is separation   angle and distance determine everything
    
    d = 2;
    th0 = 30*(pi/180);
    mm = 12;
    h = mm * (d/2); % extend cone
    a = d / ( 2 * cos(th0));  % implied a, by distance d and angle
    b = (mm*d)  / ( 2 * cos(th0)); % set it away by multiplier
    t = b * sin(th0); % get the top radius
    z0 = d/2;
    z  = +/- h/t sqrt(x^2 +y^2) +  z0         h/t = cos(th0)/sin(th0)  = cot(th0)
%}
global d; syms d; assume(d,'real'); assume(d,'positive'); % distance of needle tips
global a; syms a; assume(a, 'real'); assume(a, 'positive'); % radius of inner tip
global b; syms b; assume(b, 'real'); assume(a, 'positive'); % radius of outer tip
global th0; syms th0; assume(th0, 'real'); assume(th0,'positive'); % angle of cone
%aa = .001; % change to optimal
%bb = .002; % change to optimal
ddd = 2e-5; % change to optimal
th000 = .0558; % change to optimal
mmm = 10e3; % may have to change
aaa = ddd/(2*cos(th000));
bbb = (mmm*ddd)/(2*cos(th000));
hhh = mmm*ddd/2;

global Vs; syms Vs; assume(Vs, 'real'); % source voltage (dc) to the circuit
global Rs; syms Rs; assume(Rs, 'real'); % source resistor
global Is; syms Is; % source current
global If; syms If; % current across fuse/capacitor
global Zf; syms Zf; % impedance of the fuse/capacitor
global C; syms C; assume(C,'real'); % capacitance 
global Vf; syms Vf; % voltage across the fuse/capacitor
global Rf; syms Rf; assume(Rf, 'real'); % buffer resistor
global IL; syms IL; % current aross load
global VL; syms VL; % load voltage
global RL; syms RL; assume(RL,'real');
global I1; syms I1; % kcl/kvl holder
global I2; syms I2; % kcl/kvl holder

global s; syms s; % ilaplace() holder
global t; syms t; % laplace() holder
                
ee = cls_EE330_helper();
const = cls_CONST();



Vf_st = Vf * (log(tan(st/2)/tan((pi-th0)/2)) / log(tan(th0/2)/tan((pi-th0)/2)) );
Ef_sr_st = [ 0 , (-1/sr)*(Vf/(log(tan(th0/2)/tan((pi-th0)/2))))*(1/sin(st)) , 0 ];  
% r = a = d/2cos(th0)  th=th0 for max   ...abs()

Vmax = 12; % break down should occur at 12 V or more
Ebd = const.bd_air; % 3e6 V/m   Air  @ 1 ATM
rng_d = [.00001, .001]; % let distance "d" range from 10 um to 1 mm
rng_th0 = [pi/180, (pi/180)*89]; % let theta_0 range from 1° to 89°
buf_d = 10*rng_d(1);
buf_th0 = 10*rng_th0(1);

rng_rs = [1, 10e3];
buf_rs = 10;
rng_rf = [1, 10e3];
buf_rf = 10;
rng_rl = [1, 10e3];
buf_rl = 10;

%publisher('ee330_project.m');
%------------------------------------------------------------------------------------------ #1
if sel == 1
    fprintf('zf = oo  because it is open, steady state\n');
    I1 = Vs/(Rs+Rf+RL);
    fprintf('there is only one current, I =\n');
    pretty(I1);
    Vf = I1 * (Rf+RL);
    fprintf('voltage across fuse, Vf = \n');
    pretty(Vf);
    VL = I1 * RL;
    fprintf('voltage across load, VL = \n');
    pretty(VL);
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    fprintf('zf = 0, it is a conducting wire, the fuse was activated\n');
    fprintf('the I2 current = 0 because all current goes through the short\n');
    fprintf('Vf = 0 because everything is disipated through the Rs resistor\n');
    fprintf('the load should be protected\n');
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    v = [Vs;
          0];
    coef = [ Rs + Zf, -Zf;
            -Zf, Zf+Rf+RL];
    i = coef\v;
    I1 = i(1,1);
    I2 = i(2,1);
    simplify(I1-(I1-I2)-I2); % 0 means good
    fprintf('by KVL, I1 =\n');
    pretty(I1);
    fprintf('I2 = \n');
    pretty(I2);
    
    fprintf('\nVf = Vs - I1 Rs = I2( Rf + RL ) :\n');
    Vf = Vs - I1*Rs;
    pretty(simplify(Vf));
    fprintf('confirm Vf :\n');
    Vf = I2*(Rf+RL);
    pretty(Vf);
    
    fprintf('\nVL = I2 * RL = Vf * (RL/(RL+Rf)) = Vf - I2*Rf :\n');
    VL = I2 * RL;
    pretty(VL);
    fprintf('confirm VL:\n');
    VL = Vf - I2*Rf;
    pretty(simplify(VL));
    fprintf('confirm VL again:\n');
    VL = Vf*(RL/(RL+Rf));
    pretty(VL);
    
    Zfs = 1/(s*C); % laplace transform of fuse impedance
    Vss = Vs/s; % laplace transform of source voltage  Vs(t)= V0 u(t)
    I2s = subs(I2, [Zf,Vs], [Zfs, Vss]);
    Vfs = I2s*(Rf+RL);
    fprintf('lapace(Vf) =\n');
    pretty((Vfs));
    fprintf('ilaplace(Vfs) =\n');
    Vft = ilaplace(Vfs, s, t);
    pretty((Vft));
    VLs = I2s*RL;
    fprintf('lapace(VL) =\n');
    pretty((VLs));
    fprintf('ilaplace(VLs) =\n');
    VLt = ilaplace(VLs, s, t);
    pretty((VLt));
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    
    pts = 15;
    vpts = 5;
    buf = 1;
    rngx = 5;
    rngy = 5;
    rngz = 7;
    pos = [20, 20, 700, 700];
    ttl = 'The Fuse';
    
    %
    figure('Name', ttl,...
               'Position', pos,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    view(150,10);
    xlabel('x axis');
    ylabel('y axis'); 
    zlabel('z axis');  
    xlim([-rngx-buf, rngx+buf]);
    ylim([-rngy-buf, rngy+buf]);
    zlim([-rngz-buf, rngz+buf]);
    xax = linspace(-rngx-buf, rngx+buf , pts);
    yax = linspace(-rngy-buf, rngy+buf , pts);
    zax = linspace(-rngz-buf, rngz+buf , pts);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rngx+buf, 0          , 0           ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0          , rngy+buf, 0           ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0          , 0          , rngz+buf ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rngx+buf, 0      , 0            , '+ X');
    text(0      , rngy+buf, 0            , '+ Y');
    text(0      , 0      , rngz+buf, '+ Z');
    
    d = 1;
    th0 = 30*(pi/180);
    mm = 12;
    h = mm * (d/2); % extend cone
    a = d / ( 2 * cos(th0));  % implied a, by distance d and angle
    b = (mm*d)  / ( 2 * cos(th0)); % set it away by multiplier
    %t = b * sin(th0); % get the top radius
    z0 = d/2;
    g = cos(th0)/sin(th0);
    
    Vs = 100;
    Esph = [0,(-1/(sr*sin(st)))*(Vs/(log((tan(th0/2)/(tan((pi-th0)/2)))))),0];
    Erec = ee.transVecSR(Esph);
    
    r = linspace(a,b,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi) ./ g;
    y = R .* sin(Fi) ./ g;
    z = R;
    surf(x,y,z,'FaceColor', 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
    surf(x,y,-z,'FaceColor', 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
    plot3(0,0,d/2,'r.', 'MarkerSize', 20);
    plot3(0,0,-d/2,'r.', 'MarkerSize', 20);
    %plot3([0,0],[0,0],[d/2,-d/2], 'g-','LineWidth',3);
    
    maxE = quiver3(0,0,d/2,0, 0, -d);
    maxE.AutoScale = 'on';
    maxE.AutoScaleFactor = 0.7;
    maxE.Color = 'g';
    maxE.LineWidth = 5;
    maxE.LineStyle = '--';
    %R.Marker = '.';
    maxE.ShowArrowHead = 'on';
    
    c1 = zeros(pts,pts,pts);
    c2 = zeros(pts, pts, pts);
    c3 = zeros(pts, pts, pts);
    temp_sph = [0,0,0];
    temp_rec = [0,0,0];
    for row = 1:pts
        for col = 1:pts
            temp_rec = [x(row,col),y(row,col),z(row,col)];
            ee.feed([x(row,col),y(row,col),z(row,col)], 'R');
            temp_sph(1) = ee.pts(3,1);
            temp_sph(2) = ee.pts(3,2);
            temp_sph(3) = ee.pts(3,3);
            c1(row,col,1) = temp_rec(1);
            c2(row,col,1) = temp_rec(2);
            c3(row,col,1) = temp_rec(3);
            srn = linspace(temp_sph(1),b-intv, pts);
            stn = linspace(temp_sph(2),pi-th0-intv,pts);
            sfn = linspace(temp_sph(3), 2*pi, pts);
            for dep = 2:pts
                tsr = subs(Esph(1), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                tst = subs(Esph(2), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                tsf = subs(Esph(3), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                ee.feed([tsr, tst, tsf],'S');
               c1(row,col,dep) = ee.pts(1,1);
               c2(row,col,dep) = ee.pts(1,2);
               c3(row,col,dep) = ee.pts(1,3);
            end
        end
    end
    xfl = zeros(pts);
    yfl = zeros(pts);
    zfl = zeros(pts);
    for row = 1:pts
        for col = 1:pts
            for dep = 1:pts
                xfl(dep) = c1(row, col,dep);
                yfl(dep) = c2(row, col,dep);
                zfl(dep) = c3(row, col,dep);
                plot3(xfl(dep), yfl(dep), zfl(dep),'g.','MarkerSize', 10);
            end
            %plot3(xfl, yfl,zfl , 'g-', 'LineWidth', 1);
        end
    end
    %}
    
    %
    p = zeros(pts,pts);
    q = zeros(pts,pts);
    r = zeros(pts, pts);
    for rows = 1:pts
        for cols = 1:pts
            p(rows,cols)=subs(Erec(1), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
            q(rows,cols)=subs(Erec(2), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
            r(rows,cols)=subs(Erec(3), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
        end
    end
    
    Ef = quiver3(x,y,z,p, q, r);
    Ef.AutoScale = 'on';
    Ef.AutoScaleFactor = 2;
    Ef.Color = 'g';
    Ef.LineWidth = 3;
    Ef.LineStyle = '-';
    %R.Marker = '.';
    Ef.ShowArrowHead = 'on';
    
    
    x = linspace(-5,5,pts);
    y = linspace(-5,5,pts);
    [X,Y]=meshgrid(x,y);
    Z1 = (h/t) .* ( X.^2 + Y.^2).^(1/2) + z0;
    Z2 = (-h/t) .* ( X.^2 + Y.^2).^(1/2) - z0;
    surf(X,Y,Z1);
    surf(X,Y,Z2);
    plot3(0,0,d/2,'rx', 'MarkerSize', 10);
    plot3(0,0,-d/2,'rx', 'MarkerSize', 10);
    
    a=1;
    b=5;
    r = linspace(a,b,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi);
    y = R .* sin(Fi);
    z = R;
    surf(x,y,z);
    
    a=1;
    b=5;
    r = linspace(b,a,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi);
    y = R .* sin(Fi);
    z = -R;
    %}
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    
    pts = 35;
    pos = [20, 20, 700, 700];
    ttl = 'optimize angle and distance';
    
    dd = linspace(rng_d(1),rng_d(2), pts);
    th00 = linspace(rng_th0(1), rng_th0(2), pts);
    [D,Th] = meshgrid(dd,th00);
    V = zeros(pts,pts);
    zp = zeros(pts,pts)+Vmax;
    for row = 1:pts
        for col = 1:pts
            V(row,col) = D(row,col)*Ebd*tan(Th(row,col))*abs(log(tan(Th(row,col)/2)));
        end 
    end
    rng_v = [0,0];
    rng_v(1) = min(V, [], 'all');
    rng_v(2) = max(V, [], 'all');
    buf_v = 10 * rng_v(1);
    
    xax = linspace(rng_d(1)-buf_d    , rng_d(2)+buf_d , pts);
    yax = linspace(rng_th0(1)-buf_th0, rng_th0(2)+buf_th0 , pts);
    zax = linspace(rng_v(1)-buf_v    , rng_v(2)+buf_v , pts);
            
    %
    figure('Name', ttl,...
               'Position', pos,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    view(150,10);
    xlabel('distance (m)');
    ylabel('angle (rad)'); 
    zlabel('Vmax (V)');  
    xlim([rng_d(1)-buf_d    , rng_d(2)+buf_d]);
    ylim([rng_th0(1)-buf_th0, rng_th0(2)+buf_th0]);
    zlim([rng_v(1)-buf_v    , rng_v(2)+buf_v]);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rng_d(2)+buf_d, 0                 , 0              ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0             , rng_th0(2)+buf_th0, 0              ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0             , 0                 , rng_v(2)+buf_v ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_d(2)+buf_d, 0                 , 0             , '1 mm', 'FontSize', 16);
    text(0             , rng_th0(2)+buf_th0, 0             , '\pi/2', 'FontSize', 16);
    text(0             , 0                 , rng_v(2)+buf_v, 'V', 'FontSize', 16);
    plot3(0,0,1000, 'g.', 'MarkerSize', 20);
    text(0,0,1000, '1 kV', 'FontSize', 12);
    
    surf(D, Th, V,...
    'FaceColor', 'b',...
    'FaceAlpha', .6,...
    'EdgeColor', 'none'); 
    mesh(D, Th, V,...
    'FaceColor', 'none',...
    'EdgeColor', 'r',...
    'FaceAlpha', .9,...
    'LineStyle',':',...
    'LineWidth',.5);
    surf(D, Th, zp,...
    'FaceColor', 'k',...
    'FaceAlpha', .2,...
    'EdgeColor', 'none');
    hold off;

    figure('Name', sprintf('contours on Vmax'),...
           'Position', pos,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    %axis equal;
    view(2); % above
    title('possible solutions of angle \theta_0 and distance "d"','FontSize', 16);
    xlabel('distance');
    ylabel('angle'); 
    xlim([rng_d(1)-buf_d, rng_d(2)+buf_d]);
    ylim([rng_th0(1)-buf_th0, rng_th0(2)+buf_th0]);
    plot(xax  , 0*xax, 'k', 'linewidth', 1);
    plot(0*yax, yax  , 'k', 'linewidth', 1);
    plot(rng_d(2)+buf_d, 0                 , 'y.', 'markersize', 20, 'linewidth', 10);  
    plot(0             , rng_th0(2)+buf_th0, 'y.', 'markersize', 20, 'linewidth', 10); 
    [con_d, con_th] = contour(D, Th, V, [12,100,1000], 'g', 'LineWidth', 2);
    clabel(con_d, con_th);
    hold off;
    
    evz = 500; % adjust
    dd = linspace(rng_d(1), rng_d(2), evz);
    thh = linspace(rng_th0(1), rng_th0(2), evz);
    [DD, TH] = meshgrid(dd, thh);
    VM = zeros(evz,evz);
    theo = zeros(evz, 3);
    for row = 1:evz
        for col = 1:evz
            VM(row,col)=DD(row,col)*Ebd*tan(TH(row,col))*abs(log(tan(TH(row,col)/2))); 
        end
    end
    thtas = zeros(evz,evz);
    dists = zeros(evz,evz);
    vms = zeros(evz, evz);
    for row = 1:evz
        for col = 1:evz
            if VM(row,col)<=13
                if VM(row,col)>=11
                    thtas(row,col) = TH(row,col);
                    dists(row,col) = DD(row,col);
                    vms(row,col)=VM(row,col);
                end
            end  
        end
    end
    rng_v = [0,20];
    buf_v = 1;
    rng_d(1) = min(dists, [], 'all');
    rng_d(2) = max(dists, [], 'all');
    buf_d = 10 * rng_d(1);
    rng_th0(1) = min(thtas, [], 'all');
    rng_th0(2) = max(thtas, [], 'all');
    buf_th0 = 10 * rng_th0(1);
    xax = linspace(rng_d(1)-buf_d    , rng_d(2)+buf_d , pts);
    yax = linspace(rng_th0(1)-buf_th0, rng_th0(2)+buf_th0 , pts);
    zax = linspace(rng_v(1)-buf_v    , rng_v(2)+buf_v , pts);
    figure('Name', 'solutions',...
               'Position', pos,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    view(150,10);
    xlabel('distance (m)');
    ylabel('angle (rad)'); 
    zlabel('Vmax (V)');  
    xlim([rng_d(1)-buf_d    , rng_d(2)+buf_d]);
    ylim([rng_th0(1)-buf_th0, rng_th0(2)+buf_th0]);
    zlim([rng_v(1)-buf_v    , rng_v(2)+buf_v]);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rng_d(2)+buf_d, 0                 , 0              ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0             , rng_th0(2)+buf_th0, 0              ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0             , 0                 , rng_v(2)+buf_v ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_d(2)+buf_d, 0                 , 0             , 'max d', 'FontSize', 16);
    text(0             , rng_th0(2)+buf_th0, 0             , 'max \theta_0', 'FontSize', 16);
    text(0             , 0                 , rng_v(2)+buf_v, 'V', 'FontSize', 16);
    plot3(0,0,Vmax, 'g.', 'MarkerSize', 20);
    text(0,0,Vmax, 'Vmax = 12 V', 'FontSize', 12);
    theo_th = [999,-999];
    theo_d = [999,-999];
    for row = 1:evz
        for col = 1:evz
            if vms(row,col) ~= 0
                plot3(dists(row,col),thtas(row,col),vms(row,col),'r.','MarkerSize',15);
                if thtas(row,col)<theo_th(1) 
                    theo_th(1)=thtas(row,col);
                end
                if thtas(row,col)>theo_th(2)
                    theo_th(2)=thtas(row,col);
                end
                if dists(row,col)<theo_d(1)
                    theo_d(1)=dists(row,col);
                end
                if dists(row,col)>theo_d(2)
                    theo_d(2)=dists(row,col);
                end
            end
        end
    end
    x=linspace(theo_d(1),theo_d(2),pts);
    y=linspace(theo_th(1),theo_th(2),pts);
    [X,Y]=meshgrid(x,y);
    Z = zeros(pts,pts);
    for row = 1:pts
        for col = 1:pts
            Z(row,col)=X(row,col)*Ebd*tan(Y(row,col))*abs(log(tan(Y(row,col)/2)));
        end
    end
    surf(X, Y, Z,...
    'FaceColor', 'c',...
    'FaceAlpha', .2,...
    'EdgeColor', 'none'); 
    %mesh(X, Y, Z,...
    %'FaceColor', 'none',...
    %'EdgeColor', 'y',...
    %'FaceAlpha', .9,...
    %'LineStyle',':',...
    %'LineWidth',1);
    plot3(ddd,th000,Vmax,'go','MarkerSize',20,'LineWidth',4);
    text(ddd,th000,Vmax-2,'OPTIMAL','FontSize',20);
    fprintf('for a Vmax = %.1f :\n', Vmax);
    fprintf('theta_0 can range from %.2f to % .2f [radians]\n', theo_th(1), theo_th(2));
    fprintf(' ... %.2f to %.2f  [degrees]\n', (180/pi)*theo_th(1),(180/pi)*theo_th(2));
    fprintf('d can range from %d  to %d  [meters]\n', theo_d(1), theo_d(2));
    %}
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    t1 = (2 * pi * const.ep0);
    t2 = abs(log( (tan(th000/2)) / (tan((pi-th000)/2)) ));
    t3 = (mmm*ddd/(2*cos(th000)))-(ddd/(2*cos(th000)));
    c = (t1/t2)*t3;
    fprintf('with th0 = %.2f °  , d = %d m  , m = %.1f ;  C = %d  F\n',...
        (180/pi)*th000, ddd, mmm, double(c));
    fprintf('\n\t\tthis reults in the largest r=b= %d m\n', double((mmm*ddd)/(2*cos(th000))));
    fprintf('\n the break down voltage will occur at %.3f V\n',...
        double(ddd*Ebd*tan(th000)*abs(log(tan(th000/2)))));
    
    v_at_bd = subs(Vf_st, [Vf, st, th0], [Vmax, th000, th000]);
    ef = subs(Ef_sr_st, [sr, st, th0, Vf], [aaa, th000, th000, Vmax]);
    ef_mag = double(abs(ef(2)));
    fprintf('\nbreaking down at %.3f V, Emag = %.4f  V/m  using: th0, d, m\n', v_at_bd, ef_mag);
    fprintf('Emag must be in excess of   Ebd = %.4f  V/m\n', Ebd);
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    vs = 14;
    rs = 30;
    rf = 10;
    rl = 100;
    
    is = vs/(rs+rf+rl); % = il
    vf = vs*((rf+rl)/(rs+rf+rl));
    vl = (vs*rl)/(rs+rf+rl);
    
    fprintf('rs = %d , rf = %d , rl = %d\n', rs, rf, rl);
    fprintf('i = %d , vf = %d , vl = %d\n', is, vf, vl);
    fprintf('\nsimulate voltage surge by making rs = 0 = short\n');
    rs1 = 0;
    is = vs/(rs1+rf+rl); % = il
    vf = vs*((rf+rl)/(rs1+rf+rl));
    vl = (vs*rl)/(rs1+rf+rl);
    fprintf('\nrs = %d , rf = %d , rl = %d\n', rs1, rf, rl);
    fprintf('i = %d , vf = %d , vl = %d\n', is, vf, vl); 
    
    fprintf('\nnow the fuse activated and there is a short through the capacitor\n');
    is = vs/rs;
    fprintf('vl = vf = 0     I = %d     Vs = %d\n', is, vs);
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    pts = 15;
    vpts = 15;
    buf_m = 10;
    pos = [20, 20, 700, 700];
    ttl = 'optimized fuze';
    
    sps = [ 1*(1/5)*bbb, th000, 1*(pi/4);
            2*(1/5)*bbb, th000, 1*(pi/4);
            3*(1/5)*bbb, th000, 1*(pi/4);
            4*(1/5)*bbb, th000, 1*(pi/4);
            bbb, th000, 1*(pi/4);
            1*(1/5)*bbb, th000, 3*(pi/4);
            2*(1/5)*bbb, th000, 3*(pi/4);
            3*(1/5)*bbb, th000, 3*(pi/4);
            4*(1/5)*bbb, th000, 3*(pi/4);
            bbb, th000, 3*(pi/4);
            1*(1/5)*bbb, th000, 5*(pi/4);
            2*(1/5)*bbb, th000, 5*(pi/4);
            3*(1/5)*bbb, th000, 5*(pi/4);
            4*(1/5)*bbb, th000, 5*(pi/4);
            bbb, th000, 5*(pi/4);
            1*(1/5)*bbb, th000, 7*(pi/4);
            2*(1/5)*bbb, th000, 7*(pi/4);
            3*(1/5)*bbb, th000, 7*(pi/4);
            4*(1/5)*bbb, th000, 7*(pi/4);
            bbb, th000, 7*(pi/4) ];
        
    eps = [ 1*(1/5)*bbb, pi-th000, 1*(pi/4);
            2*(1/5)*bbb, pi-th000, 1*(pi/4);
            3*(1/5)*bbb, pi-th000, 1*(pi/4);
            4*(1/5)*bbb, pi-th000, 1*(pi/4);
            bbb, pi-th000, 1*(pi/4);
            1*(1/5)*bbb, pi-th000, 3*(pi/4);
            2*(1/5)*bbb, pi-th000, 3*(pi/4);
            3*(1/5)*bbb, pi-th000, 3*(pi/4);
            4*(1/5)*bbb, pi-th000, 3*(pi/4);
            bbb, pi-th000, 3*(pi/4);
            1*(1/5)*bbb, pi-th000, 5*(pi/4);
            2*(1/5)*bbb, pi-th000, 5*(pi/4);
            3*(1/5)*bbb, pi-th000, 5*(pi/4);
            4*(1/5)*bbb, pi-th000, 5*(pi/4);
            bbb, pi-th000, 5*(pi/4);
            1*(1/5)*bbb, pi-th000, 7*(pi/4);
            2*(1/5)*bbb, pi-th000, 7*(pi/4);
            3*(1/5)*bbb, pi-th000, 7*(pi/4);
            4*(1/5)*bbb, pi-th000, 7*(pi/4);
            bbb, pi-th000, 7*(pi/4) ];    
    
            
    Evec = zeros(20, 3, vpts+2);
    for row = 1:20
        Evec(row,1,1) = sps(row,1);
        Evec(row,2,1) = sps(row,2);
        Evec(row,3,1) = sps(row,3);
        Evec(row,1,vpts+2) = eps(row,1);
        Evec(row,2,vpts+2) = eps(row,2);
        Evec(row,3,vpts+2) = eps(row,3);
    end
    fix = .02;
    tran_th = linspace(th000+fix,pi-th000-fix,vpts);
    for row = 1:20
        for stp = 2:vpts+1
            Evec(row,1,stp) = Evec(row,1,1);
            Evec(row,2,stp) = tran_th(stp-1);
            Evec(row,3,stp) = Evec(row,3,1);
        end
    end
    for row = 1:20
        for dep = 1:vpts+2
            tmp = [ Evec(row,1,dep), Evec(row,2,dep), Evec(row,3,dep) ];
            ee.feed(tmp, 'S');
            Evec(row,1,dep) = ee.pts(1,1);
            Evec(row,2,dep) = ee.pts(1,2);
            Evec(row,3,dep) = ee.pts(1,3);
        end
    end
    
    rs = linspace(aaa,bbb, pts);
    fis = linspace(0,2*pi,pts);
    [R, Fi] = meshgrid(rs, fis);
    X = zeros(pts,pts);
    Y = X;
    Z = X;
    for row = 1:pts
        for col = 1:pts
            ee.feed([ R(row,col), th000, Fi(row,col)], 'S');
            X(row,col) = ee.pts(1,1);
            Y(row,col) = ee.pts(1,2);
            Z(row,col) = ee.pts(1,3);
        end
    end
    rngz = [ -hhh, hhh];
    bufz = ddd*buf_m;
    rngx = rngz;
    bufx = bufz;
    rngy = rngz;
    bufy = bufz;

    %
    figure('Name', ttl,...
               'Position', pos,...
               'NumberTitle', 'off');
    hold on;
    grid on;
    axis manual;
    view(150,10);
    xlabel('x axis [m]');
    ylabel('y axis [m]'); 
    zlabel('z axis [m]');  
    xlim([double(rngx(1)-bufx), double(rngx(2)+bufx) ]);
    ylim([double(rngy(1)-bufy), double(rngy(2)+bufy) ]);
    zlim([double(rngz(1)-bufz), double(rngz(2)+bufz) ]);
    xax = linspace(rngx(1)-bufx, rngx(2)+bufx, pts);
    yax = linspace(rngy(1)-bufy, rngy(2)+bufy, pts);
    zax = linspace(rngz(1)-bufz, rngz(2)+bufz, pts);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rngx(2)+bufx, 0          , 0           ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0          , rngy(2)+bufy, 0           ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0          , 0          , rngz(2)+bufz ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rngx(2)+bufx, 0           , 0           , '+ X');
    text(0           , rngy(2)+bufy, 0           , '+ Y');
    text(0           , 0           , rngz(2)+bufz, '+ Z');
    
    Erec = ee.transVecSR(Ef_sr_st);
    for row = 1:20
        init = [ Evec(row,1,1), Evec(row,2,1), Evec(row,3,1) ];
        eval_e = subs(Erec, [rx,ry,rz,Vf,th0], [init(1),init(2),init(3),Vmax,th000]);
        quiver3(init(1), init(2), init(3), eval_e(1), eval_e(2), eval_e(3),...
            'AutoScale', 'on',...
            'AutoScaleFactor', .00005,...
            'Color', 'g',...
            'LineWidth', 2,...
            'LineStyle', '-',...
            'ShowArrowHead', 'on');
        for ef = 1:vpts+1
            pre = [ Evec(row,1,ef), Evec(row,2,ef), Evec(row,3,ef) ];
            nxt = [ Evec(row,1,ef+1), Evec(row,2,ef+1), Evec(row,3,ef+1) ];
            stp = nxt-pre;
            quiver3(pre(1), pre(2), pre(3), stp(1), stp(2), stp(3),...
                'AutoScale', 'on',...
                'AutoScaleFactor', .7,...
                'Color', 'r',...
                'LineWidth', 2,...
                'LineStyle', '--',...
                'ShowArrowHead', 'on');
        end
    end
                   
    surf(X,Y,Z,'FaceColor', 'k', 'FaceAlpha', .8, 'EdgeColor', 'none');
    plot3(0,0,ddd/2, 'b.', 'MarkerSize', 15);
    surf(X,Y,-Z,'FaceColor', 'k', 'FaceAlpha', .8, 'EdgeColor', 'none');
    plot3(0,0,-ddd/2, 'b.', 'MarkerSize', 10);
    hold off;
    %{
    
    mesh(X, Y, Z,...
        'FaceColor', 'k',...
        'EdgeColor', 'k',...
        'FaceAlpha', .1,...
        'LineStyle',':',...
        'LineWidth',1);
    
    x = linspace(-bbb*sin(th000), bbb*sin(th000), pts);
    y = x;
    [X,Y] = meshgrid(x,y);
    Zt = zeros(pts, pts);
    Zb = Zt;
    for row = 1:pts
        for col = 1:pts
            Zt(row,col) = cot(th000)*sqrt((X(row,col)^2 + Y(row,col)^2)) + (ddd/2);
            Zb(row,col) = (-1)*cot(th000)*sqrt((X(row,col)^2 + Y(row,col)^2)) - (ddd/2);
        end
    end
    rngx = [ -bbb*sin(th000), bbb*sin(th000) ];
    bufx = aaa*buf_m;
    rngy = rngx;
    bufy = bufx;
    rngz = [ -hhh, hhh];
    bufz = ddd*buf_m;
    
    x = linspace(aaa*sin(th000), bbb*sin(th000), pts);
    y = x;
    [X,Y] = meshgrid(x,y);
    Zt = zeros(pts, pts);
    Zb = Zt;
    for row = 1:pts
        for col = 1:pts
            Zt(row,col) = cot(th000)*sqrt((X(row,col)^2 + Y(row,col)^2)) + (ddd/2);
            Zb(row,col) = (-1)*cot(th000)*sqrt((X(row,col)^2 + Y(row,col)^2)) - (ddd/2);
        end
    end
    rngx = [ -bbb*sin(th000), bbb*sin(th000) ];
    bufx = aaa*buf_m;
    rngy = rngx;
    bufy = bufx;
    rngz = [ -hhh, hhh];
    bufz = ddd*buf_m;
    
    d = 1;
    th0 = 30*(pi/180);
    mm = 12;
    h = mm * (d/2); % extend cone
    a = d / ( 2 * cos(th0));  % implied a, by distance d and angle
    b = (mm*d)  / ( 2 * cos(th0)); % set it away by multiplier
    %t = b * sin(th0); % get the top radius
    z0 = d/2;
    g = cos(th0)/sin(th0);
    
    Vs = 100;
    Esph = [0,(-1/(sr*sin(st)))*(Vs/(log((tan(th0/2)/(tan((pi-th0)/2)))))),0];
    Erec = ee.transVecSR(Esph);
    
    r = linspace(a,b,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi) ./ g;
    y = R .* sin(Fi) ./ g;
    z = R;
    surf(x,y,z,'FaceColor', 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
    surf(x,y,-z,'FaceColor', 'k', 'FaceAlpha', .2, 'EdgeColor', 'none');
    plot3(0,0,d/2,'r.', 'MarkerSize', 20);
    plot3(0,0,-d/2,'r.', 'MarkerSize', 20);
    %plot3([0,0],[0,0],[d/2,-d/2], 'g-','LineWidth',3);
    
    maxE = quiver3(0,0,d/2,0, 0, -d);
    maxE.AutoScale = 'on';
    maxE.AutoScaleFactor = .7;
    maxE.Color = 'g';
    maxE.LineWidth = 5;
    maxE.LineStyle = '--';
    %R.Marker = '.';
    maxE.ShowArrowHead = 'on';
    
    c1 = zeros(pts,pts,pts);
    c2 = zeros(pts, pts, pts);
    c3 = zeros(pts, pts, pts);
    temp_sph = [0,0,0];
    temp_rec = [0,0,0];
    for row = 1:pts
        for col = 1:pts
            temp_rec = [x(row,col),y(row,col),z(row,col)];
            ee.feed([x(row,col),y(row,col),z(row,col)], 'R');
            temp_sph(1) = ee.pts(3,1);
            temp_sph(2) = ee.pts(3,2);
            temp_sph(3) = ee.pts(3,3);
            c1(row,col,1) = temp_rec(1);
            c2(row,col,1) = temp_rec(2);
            c3(row,col,1) = temp_rec(3);
            srn = linspace(temp_sph(1),b-intv, pts);
            stn = linspace(temp_sph(2),pi-th0-intv,pts);
            sfn = linspace(temp_sph(3), 2*pi, pts);
            for dep = 2:pts
                tsr = subs(Esph(1), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                tst = subs(Esph(2), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                tsf = subs(Esph(3), [sr, st, sf], [srn(dep), stn(dep), sfn(dep)]);
                ee.feed([tsr, tst, tsf],'S');
               c1(row,col,dep) = ee.pts(1,1);
               c2(row,col,dep) = ee.pts(1,2);
               c3(row,col,dep) = ee.pts(1,3);
            end
        end
    end
    xfl = zeros(pts);
    yfl = zeros(pts);
    zfl = zeros(pts);
    for row = 1:pts
        for col = 1:pts
            for dep = 1:pts
                xfl(dep) = c1(row, col,dep);
                yfl(dep) = c2(row, col,dep);
                zfl(dep) = c3(row, col,dep);
                plot3(xfl(dep), yfl(dep), zfl(dep),'g.','MarkerSize', 10);
            end
            %plot3(xfl, yfl,zfl , 'g-', 'LineWidth', 1);
        end
    end
    
    p = zeros(pts,pts);
    q = zeros(pts,pts);
    r = zeros(pts, pts);
    for rows = 1:pts
        for cols = 1:pts
            p(rows,cols)=subs(Erec(1), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
            q(rows,cols)=subs(Erec(2), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
            r(rows,cols)=subs(Erec(3), [rx, ry, rz], [x(rows,cols),y(rows,cols),z(rows,cols)]);
        end
    end
    
    Ef = quiver3(x,y,z,p, q, r);
    Ef.AutoScale = 'on';
    Ef.AutoScaleFactor = 2;
    Ef.Color = 'g';
    Ef.LineWidth = 3;
    Ef.LineStyle = '-';
    %R.Marker = '.';
    Ef.ShowArrowHead = 'on';
    
    
    x = linspace(-5,5,pts);
    y = linspace(-5,5,pts);
    [X,Y]=meshgrid(x,y);
    Z1 = (h/t) .* ( X.^2 + Y.^2).^(1/2) + z0;
    Z2 = (-h/t) .* ( X.^2 + Y.^2).^(1/2) - z0;
    surf(X,Y,Z1);
    surf(X,Y,Z2);
    plot3(0,0,d/2,'rx', 'MarkerSize', 10);
    plot3(0,0,-d/2,'rx', 'MarkerSize', 10);
    
    a=1;
    b=5;
    r = linspace(a,b,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi);
    y = R .* sin(Fi);
    z = R;
    surf(x,y,z);
    
    a=1;
    b=5;
    r = linspace(b,a,pts);
    fi = linspace(0,2*pi,pts);
    [R,Fi]=meshgrid(r,fi);
    x = R .* cos(Fi);
    y = R .* sin(Fi);
    z = -R;
    %}
end
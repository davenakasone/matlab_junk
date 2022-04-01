%{
    ee360 ch10

    #0  get user-defined functions correct

    #1021.3  hw7, 10.21c
    #1021.4  hw7, 10.21d
    #1021.5  hw7, 10.21e
    #1021.6  hw7, 10.21f
    #1021.7  hw7, 10.21g
    #1028    hw7, 10.28
    #1034    hw7, 10.34
    #1037    hw7, 10.37
    #1042.1  hw7, 1042.1

    #772 chw5, #2 magnitude ->  Mesh C
    #773 chw5, #3 freq resp (mag and phase)
    #774 chw5, #4 roots, poles, zeros, and system function  + #5 pol/zero
    #776 chw5, #6 pfe, plot h[n]  0<=n<=50
    #777   p7


%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
old_val = sympref('HeavisideAtOrigin', 1);
format shortE; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
%format compact; % [compact,loose]
                
global Z; syms Z;            % your temporary Z for integration and differentiation  * also Z trans
global null; null = double.empty();  % nulled array of size 1x1
global intv; intv = eps*(1e14);
global in; syms in; assume(in, {'real', 'integer'}); % index n

ee = cls_EE330_helper();
const = cls_CONST();
%publisher('ee330_hw4.m');


                            sel = 777;  % CHANGE CHANGE CHANGE


%------------------------------------------------------------------------------------------ #0
if sel == 0
    hz = ( (1+(1/(3*Z))) ) / ( (1-(1/(2*Z)))    );
    %pretty(simplify(hz));
    tst = fun_pfe_ztran(hz);
    fun_DT_zer_pol(zers, pols, roc, side);
end


%------------------------------------------------------------------------------------------ #1021.3
if sel == 1021.3
    xz = 1 /  (1-(1/Z))    ;
    zp = fun_pfe_ztran(xz);
    fun_DT_zer_pol(zp(1), zp(2), 1, 999);
end


%------------------------------------------------------------------------------------------ #1021.4
if sel == 1021.4
    xz = 4*Z^3 * (1 /  (1-(1/(2*Z))))    ;
    zp = fun_pfe_ztran(xz);
    %fun_DT_zer_pol(zp(1), zp(2), 1, 999);
    fun_DT_zer_pol(0, .5, 1, 999);
end


%------------------------------------------------------------------------------------------ #1021.5
if sel == 1021.5
    xz = 9*Z^2 * (1 /  (1+3*Z) )    ;
    zp = fun_pfe_ztran(xz);
    %fun_DT_zer_pol(zp(1), zp(2), 1, 999);
    fun_DT_zer_pol([0,0], -1/3, 2, 999);
end


%------------------------------------------------------------------------------------------ #1021.6
if sel == 1021.6
    xz = (1/(4*Z))^3 * (1/(4*Z-1));
    zp = fun_pfe_ztran(xz);
    %fun_DT_zer_pol(zp(1), zp(2), 2, 999);
    fun_DT_zer_pol(null, [0,1/4], 3, [.05,1/4]);
end


%------------------------------------------------------------------------------------------ #1021.7
if sel == 1021.7
    x1 = (2/(2-Z));
    x2 = (1/(4*Z-1));
    xz = x1*x2;
    zp = fun_pfe_ztran(xz);
    %fun_DT_zer_pol(zp(1), zp(2), 2, 999);
    fun_DT_zer_pol(null, [1/4,2], 3, [1/4,2]);
end


%------------------------------------------------------------------------------------------ #1028
if sel == 1028
    xz = 1-(.95/Z^6);
    [num,den] = numden(xz);
    zer = coeffs(num, Z, 'all');
    pol = coeffs(den, Z, 'all');
    rt_num = double(transpose(roots(zer)));
    rt_den = double(transpose(roots(pol)));
    %p = fun_pfe_ztran(xz);
    %fun_DT_zer_pol(zp(1), zp(2), 2, 999);
    %fun_DT_zer_pol(rt_num, rt_den, 999, [.2,2]);
    
    pos = [20,20,700,700];
    pts = 200;
    rng_w = [0,2*pi];
    buf_w = 1;
    rng_mag = [0,0];
    buf_m = 1;
    rng_phs = [0,0];
    buf_p = 1;
    syms s;
    w = linspace(rng_w(1),rng_w(2),pts);
    Hss = subs(xz, Z, s);
    [num,den] = numden(Hss);
    b_cof = double(coeffs(num,s,'All'));
    a_cof = double(coeffs(den,s,'All'));
    H = freqs(b_cof, a_cof, w);
    
    Hmag = zeros(1,pts);
    Hphs = zeros(1,pts);
    for ct = 1:pts
        Hmag(1,ct) = 20*log10(abs(H(ct)));
        Hphs(1,ct) = (180/pi)*angle(H(ct));
        if Hmag(1,ct) < rng_mag(1) 
            rng_mag(1) = Hmag(1,ct);
        end
        if Hmag(1,ct) > rng_mag(2)
            rng_mag(2) = Hmag(1,ct);
        end
        if Hphs(1,ct) < rng_phs(1) 
            rng_phs(1) = Hphs(1,ct);
        end
        if Hphs(1,ct) > rng_phs(2)
            rng_phs(2) = Hphs(1,ct);
        end
    end
    ax_w = linspace(rng_w(1)-buf_w, rng_w(2)+buf_w, pts);
    ax_mag = linspace(rng_mag(1)-buf_m, rng_mag(2)+buf_m, pts);
    ax_phs = linspace(rng_phs(1)-buf_p, rng_phs(2)+buf_p, pts);
    figure('Position', pos);
    
    subplot(2,1,1);
    hold on;
    grid on;
    title('MAGNITUDE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('H dB = 20log10(|H|');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_mag(1)-buf_m, rng_mag(2)+buf_m]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_mag, ax_mag  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_mag(2)+buf_m, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hmag,'r-', 'LineWidth', 2);
    
    subplot(2,1,2);
    hold on;
    grid on;
    title('PHASE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('phase in ¬∞');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_phs(1)-buf_p, rng_phs(2)+buf_p]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_phs, ax_phs  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_phs(2)+buf_p, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hphs,'r-', 'LineWidth', 2);
end

%------------------------------------------------------------------------------------------ #1034
if sel == 1034
    h1 = 1/Z;
    h2 = 1-(1/Z)-(1/Z^2);
    hz = h1/h2;
    %pretty(simplify(hz));
    %tst = fun_pfe_ztran(hz);
    [num,den]= numden(hz);
    num_cof = coeffs(num,Z,'all');
    den_cof = coeffs(den, Z, 'all');
    num_rot = double(transpose(roots(num_cof)));
    den_rot = double(transpose(roots(den_cof)));
    fun_DT_zer_pol(num_rot, den_rot, 3,[den_rot(1),den_rot(2)]);
end


%------------------------------------------------------------------------------------------ #1037
if sel == 1037
    h1 = 1-(9/(8*Z));
    h2 = 1+(1/(3*Z))-(2/(9*Z^2));
    hz = h1/h2;
    pretty(simplify(hz));
    %tst = fun_pfe_ztran(hz);
    [num,den]= numden(hz);
    num_cof = coeffs(num,Z,'all');
    den_cof = coeffs(den, Z, 'all');
    num_rot = double(transpose(roots(num_cof)));
    den_rot = double(transpose(roots(den_cof)));
    fun_DT_zer_pol(num_rot, den_rot, 1,[den_rot(1),den_rot(2)]);
end


%------------------------------------------------------------------------------------------ #1042.1
if sel == 1042.1
    h1 = 1-(1/Z);
    h2 = 1-(3/(4*Z))+(1/(8*Z^2))    ;
    hz = h1/h2;
    %pretty(simplify(hz));
    tst = fun_pfe_ztran(hz);
    [num,den]= numden(hz);
    num_cof = coeffs(num,Z,'all');
    den_cof = coeffs(den, Z, 'all');
    num_rot = double(transpose(roots(num_cof)));
    den_rot = double(transpose(roots(den_cof)));
    %fun_DT_zer_pol(num_rot, den_rot, 1,[den_rot(1),den_rot(2)]);
end


r = .9;
w0 = pi/7;
Hz = ( 1 - ((r * cos(w0))/Z) ) / ( 1 - ( ( 2 * r * cos(w0) ) / Z) + (r^2 / Z^2) );
Hz = expand(Hz);
[h_num, h_den] = numden(Hz);
h_num_cof = coeffs(h_num, Z, 'all');
h_den_cof = coeffs(h_den, Z, 'all');
h_num_rot = double(transpose(roots(h_num_cof)));
h_den_rot = double(transpose(roots(h_den_cof)));
%------------------------------------------------------------------------------------------ #772
if sel == 772
    
    bufS = 1;
    bufO = 1;
    bufH = 1;
    pts = 25;
    fil = 128;
    pos = [20,20,700,700];
    rng_sig = [-2,2];
    rng_omg = [-2,2];
    rng_Hmag = [0,0];
    tis = '20 log10 ( Hmag ) vs z real[-2,2] imag[-2,2]';
    xax = linspace(rng_sig(1)-bufS, rng_sig(2)+bufS, fil);
    yax = linspace(rng_omg(1)-bufO, rng_omg(2)+bufO, fil);
    sig_ax = linspace(rng_sig(1)-bufS, rng_sig(2)+bufS, pts);
    omg_ax = linspace(rng_omg(1)-bufO, rng_omg(2)+bufO, pts);
    
    Hdb = zeros(pts,pts);
    [Sig,Omg] = meshgrid(sig_ax, omg_ax);
    for row = 1:pts
        for col = 1:pts
            temp_var = Sig(row,col)+ 1j* Omg(row,col);
            Hdb(row,col)= 20*log10(abs(subs(Hz, Z, temp_var+intv)));
            if Hdb(row,col)< rng_Hmag(1)
                rng_Hmag(1) = Hdb(row,col);
            end
            if Hdb(row,col) > rng_Hmag(2)
                rng_Hmag(2) = Hdb(row,col);
            end
        end
    end
    zax = linspace(rng_Hmag(1)-bufH, rng_Hmag(2)+bufH, fil);
    
    figure('Name', 'question 2',...
           'Position', pos,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    view(125,30); % CHANGE
    title(tis, 'FontSize', 16);
    xlabel('real(z)' , 'FontSize', 14);
    ylabel('imag(z)' , 'FontSize', 14); 
    zlabel('Hdb = 20log10(|H|', 'FontSize', 14);
    xlim([rng_sig(1)-bufS, rng_sig(2)+bufS]);
    ylim([rng_omg(1)-bufO, rng_omg(2)+bufO]);
    zlim([rng_Hmag(1)-bufH, rng_Hmag(2)+bufH]);
    plot3(xax  , 0*xax, 0*xax, 'k', 'linewidth', 1);
    plot3(0*yax, yax  , 0*yax, 'k', 'linewidth', 1);
    plot3(0*zax, 0*zax, zax  , 'k', 'linewidth', 1);
    plot3(rng_sig(2)+bufS , 0          , 0                ,'y.', 'markersize', 20, 'linewidth', 10);  
    plot3(0           , rng_omg(2)+bufO, 0                ,'y.', 'markersize', 20, 'linewidth', 10); 
    plot3(0           , 0              , rng_Hmag(2)+bufH ,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_sig(2)+bufS, 0              , 0               , '+ real(z)', 'FontSize', 14);
    text(0              , rng_omg(2)+bufO, 0               , '+ imag(z)', 'FontSize', 14);
    text(0              , 0              , rng_Hmag(2)+bufH, '+ |H| dB', 'FontSize', 14);
    
    for ct = 1:length(h_num_rot)
        plot3(real(h_num_rot(ct)), imag(h_num_rot(ct)), rng_Hmag(1), 'bo', 'LineWidth', 3);
        zstr = sprintf('%.2f , %.2f', real(h_num_rot(ct)), imag(h_num_rot(ct)));
        text(real(h_num_rot(ct))+.1, imag(h_num_rot(ct))+.1, rng_Hmag(1), zstr, 'FontSize', 10);
        for ind = 1:fil
            plot3(real(h_num_rot(ct)), imag(h_num_rot(ct)), zax(ct), 'k.', 'MarkerSize', 1);
        end
    end
    
    for ct = 1:length(h_den_rot)
        plot3(real(h_den_rot(ct)), imag(h_den_rot(ct)), rng_Hmag(2), 'rx', 'LineWidth', 3);
        pstr = sprintf('%.2f , %.2f', real(h_den_rot(ct)), imag(h_den_rot(ct)));
        text(real(h_den_rot(ct))+.1, imag(h_den_rot(ct))+.1, rng_Hmag(2), pstr, 'FontSize', 10);
        for ind = 1:fil
            plot3(real(h_den_rot(ct)), imag(h_den_rot(ct)), zax(ct), 'k.', 'MarkerSize', 1);
        end
    end
    
    colormap jet;
    colorbar;
    hand = meshc(Sig,Omg,Hdb, 'FaceAlpha', .3, 'LineWidth', 1, 'FaceColor', 'texturemap');
    hand(2).LineWidth = 2;
end


%------------------------------------------------------------------------------------------ #773
if sel == 773
    pos = [20,20,700,700];
    pts = 200;
    rng_w = [-2*pi,2*pi];
    buf_w = 1;
    rng_mag = [0,0];
    buf_m = 1;
    rng_phs = [0,0];
    buf_p = 1;
    
    w = linspace(rng_w(1),rng_w(2),pts);
    H = freqs(double(h_num_cof), double(h_den_cof), w);
    
    Hmag = zeros(1,pts);
    Hphs = zeros(1,pts);
    for ct = 1:pts
        Hmag(1,ct) = 20*log10(abs(H(ct)));
        Hphs(1,ct) = (180/pi)*angle(H(ct));
        if Hmag(1,ct) < rng_mag(1) 
            rng_mag(1) = Hmag(1,ct);
        end
        if Hmag(1,ct) > rng_mag(2)
            rng_mag(2) = Hmag(1,ct);
        end
        if Hphs(1,ct) < rng_phs(1) 
            rng_phs(1) = Hphs(1,ct);
        end
        if Hphs(1,ct) > rng_phs(2)
            rng_phs(2) = Hphs(1,ct);
        end
    end
    ax_w = linspace(rng_w(1)-buf_w, rng_w(2)+buf_w, pts);
    ax_mag = linspace(rng_mag(1)-buf_m, rng_mag(2)+buf_m, pts);
    ax_phs = linspace(rng_phs(1)-buf_p, rng_phs(2)+buf_p, pts);
    figure('Position', pos);
    
    subplot(2,1,1);
    hold on;
    grid on;
    title('MAGNITUDE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('H dB = 20log10(|H|');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_mag(1)-buf_m, rng_mag(2)+buf_m]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_mag, ax_mag  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_mag(2)+buf_m, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hmag,'r-', 'LineWidth', 2);
    
    subplot(2,1,2);
    hold on;
    grid on;
    title('PHASE');
    view(2); 
    xlabel('\omega', 'FontSize', 16);
    ylabel('phase in degrees');   
    xlim([rng_w(1)-buf_w, rng_w(2)+buf_w]);
    ylim([rng_phs(1)-buf_p, rng_phs(2)+buf_p]);
    plot(ax_w      , 0 * ax_w, 'k', 'linewidth', 1);
    plot(0 * ax_phs, ax_phs  , 'k', 'linewidth', 1);        
    plot(rng_w(2)+buf_w, 0      , 'y.', 'markersize', 20, 'linewidth', 10); % +x / inputs
    plot(0      , rng_phs(2)+buf_p, 'y.', 'markersize', 20, 'linewidth', 10); % +y / outputs
    plot(w,Hphs,'r-', 'LineWidth', 2);
end


%------------------------------------------------------------------------------------------ #774
if sel == 774
    fprintf('zero 1 of 2 : %.3f + j %.3f\n', real(h_num_rot(1)), imag(h_num_rot(1)));
    fprintf('\nzero 2 of 2 : %.3f + j %.3f\n', real(h_num_rot(2)), imag(h_num_rot(2)));
    fprintf('\npole 1 of 2 : %.3f - j %.3f\n', real(h_den_rot(1)), imag(h_den_rot(1)));
    fprintf('\npole 2 of 2 : %.3f + j %.3f\n', real(h_den_rot(2)), imag(h_den_rot(2)));
    
    fprintf('\nthe system function H(z): \n\n');
    %pretty(simplify(Hz));
    
    zplane(roots(h_num_cof), roots(h_den_cof), 'g');
    fun_DT_zer_pol(h_num_rot, h_den_rot, 1, 999);
    
    %{
    b_cof = double(coeffs(num,s,'All'));
    b_rot = roots(b_cof)
    b_pol = poly(b_rot);
    a_cof = double(coeffs(den,s,'All'));
    a_rot = roots(a_cof)
    a_pol = poly(a_rot);
    hand = zplane(b_pol, a_pol, 'g');
    
    h_num_cof = coeffs(h_num, Z, 'all');
    h_den_cof = coeffs(h_den, Z, 'all');
    h_num_rot = double(transpose(roots(h_num_cof)))
    h_den_rot = double(transpose(roots(h_den_cof)))
    %}
end


%------------------------------------------------------------------------------------------ #776
if sel == 776
    %fprintf('original H(z):\n\n');
    %pretty(Hz);
    [r,p,k] = residuez(double(h_num_cof), double(h_den_cof));
    H_pfe = r(1,1)/(1-(p(1,1)/Z)) + r(2,1)/(1-(p(2,1)/Z));
    %fprintf('pratial fraction expansion of H(z) using +/- eps accuracy:\n\n');
    %pretty(H_pfe);
    %check = simplify(subs(H_pfe-Hz, Z, 1));
    %fprintf('the difference between forms is %.3f\n', check);
    hn = iztrans(Hz,Z,in);
    hn_theo = r(1,1)*(p(1,1))^in * heaviside(in) + r(2,1)*(p(2,1))^in * heaviside(in);
    tstt = double(subs(hn-hn_theo, in, 1));
    
    %
    buf_n = 1;
    buf_rl = 1;
    buf_im = 1;
    rng_h_rl = [999,-999];
    rng_h_im = [999,-999];
    rng_n = [0,50];
    fills = 128;
    pts = rng_n(2)-rng_n(1)+1;
    pos = [20,20,700,700];
    
    h_rl = zeros(1,pts);
    h_im = zeros(1,pts);
    n_in = linspace(rng_n(1), rng_n(2), pts);
    
    for ct = 1:pts
        hn_out = subs(hn_theo, in, n_in(1,ct));
        h_rl(1,ct) = real(hn_out);
        h_im(1,ct) = imag(hn_out);
        if h_rl(1,ct) < rng_h_rl(1)
            rng_h_rl(1) = h_rl(1,ct);
        end
        if h_rl(1,ct) > rng_h_rl(2)
            rng_h_rl(2) = h_rl(1,ct);
        end
        if h_im(1,ct) < rng_h_im(1)
            rng_h_im(1) = h_im(1,ct);
        end
        if h_im(1,ct) > rng_h_im(2)
            rng_h_im(2) = h_im(1,ct);
        end
    end
    
    ax_h_rl = linspace(rng_h_rl(1)-buf_rl, rng_h_rl(2)+buf_rl, fills);
    ax_h_im = linspace(rng_h_im(1)-buf_im, rng_h_im(2)+buf_im, fills);
    ax_n = linspace(rng_n(1)-buf_n, rng_n(2)+buf_n, fills);
    
    figure('Name', 'problem6',...
           'Position', pos,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    view(125,30); % CHANGE
    
    title('h[n]  0 <= n <= 50', 'FontSize', 16);
    xlabel('input n', 'FontSize', 16);
    ylabel('imag  output h[n]', 'FontSize', 16);
    zlabel('real  output h[n]', 'FontSize', 16);
    xlim([rng_n(1)-buf_n, rng_n(2)+buf_n]);
    ylim([rng_h_im(1)-buf_im, rng_h_im(2)+buf_im]);
    zlim([rng_h_rl(1)-buf_rl, rng_h_rl(2)+buf_rl]);
    plot3(ax_n     , 0*ax_n   , 0*ax_n   , 'k', 'linewidth', 1);
    plot3(0*ax_h_im, ax_h_im  , 0*ax_h_im  , 'k', 'linewidth', 1);
    plot3(0*ax_h_rl, 0*ax_h_rl, ax_h_rl, 'k', 'linewidth', 1);
    plot3(rng_n(2)+buf_n  , 0                 , 0                 ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0               , rng_h_im(2)+buf_im, 0                 ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0               , 0                 , rng_h_rl(2)+buf_rl,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_n(2)+buf_n, 0                 , 0                 , 'n in', 'FontSize', 16);
    text(0             , rng_h_im(2)+buf_im, 0                 , 'imag out', 'FontSize', 16);
    text(0             , 0                 , rng_h_rl(2)+buf_rl, 'real out', 'FontSize', 16);
    
    stmH = stem3(n_in, h_im, h_rl);
    stmH.LineStyle = '-';
    set(stmH, 'Marker', '^', 'MarkerSize', 5);
    stmH.LineWidth = 2;
    stmH.MarkerFaceColor = 'b';
    stmH.MarkerEdgeColor = 'b';
    stmH.Color = 'r';
    temH = stmH.BaseLine;
    temH.Visible = 'off';  
    %}
end


%------------------------------------------------------------------------------------------ #777
if sel == 777
    
    xn = heaviside(in);
    Xz = ztrans(xn, in, Z);
    Yz = Xz*Hz;
    [y_num, y_den] = numden(Yz);
    y_num_cof = double(coeffs(y_num));
    y_den_cof = double(coeffs(y_den));
    %fprintf('expanded Y(z):\n\n');
    %pretty(Yz);
    [r,p,k] = residuez(y_num_cof, y_den_cof);
    Y_pfe = r(1,1)/(1-(p(1,1)/Z)) + r(2,1)/(1-(p(2,1)/Z)) + r(3,1)/(1-(p(3,1)/Z));
    %fprintf('pratial fraction expansion of Y(z) using +/- eps accuracy:\n\n');
    %pretty(Y_pfe);
    %check = simplify(subs(Y_pfe-Yz, Z, 1));
    %fprintf('the difference between forms is %.3f\n', check);
    %hn = iztrans(Hz,Z,in);
    %hn_theo = r(1,1)*(p(1,1))^in * heaviside(in) + r(2,1)*(p(2,1))^in * heaviside(in);
    %tstt = double(subs(hn-hn_theo, in, 1));
    hn_theo = iztrans(Yz, Z, in);
    %
    buf_n = 1;
    buf_rl = 1;
    buf_im = 1;
    rng_h_rl = [999,-999];
    rng_h_im = [999,-999];
    rng_n = [0,50];
    fills = 128;
    pts = rng_n(2)-rng_n(1)+1;
    pos = [20,20,700,700];
    
    h_rl = zeros(1,pts);
    h_im = zeros(1,pts);
    n_in = linspace(rng_n(1), rng_n(2), pts);
    
    for ct = 1:pts
        hn_out = subs(hn_theo, in, n_in(1,ct));
        h_rl(1,ct) = real(hn_out);
        h_im(1,ct) = imag(hn_out);
        if h_rl(1,ct) < rng_h_rl(1)
            rng_h_rl(1) = h_rl(1,ct);
        end
        if h_rl(1,ct) > rng_h_rl(2)
            rng_h_rl(2) = h_rl(1,ct);
        end
        if h_im(1,ct) < rng_h_im(1)
            rng_h_im(1) = h_im(1,ct);
        end
        if h_im(1,ct) > rng_h_im(2)
            rng_h_im(2) = h_im(1,ct);
        end
    end
    
    ax_h_rl = linspace(rng_h_rl(1)-buf_rl, rng_h_rl(2)+buf_rl, fills);
    ax_h_im = linspace(rng_h_im(1)-buf_im, rng_h_im(2)+buf_im, fills);
    ax_n = linspace(rng_n(1)-buf_n, rng_n(2)+buf_n, fills);
    
    figure('Name', 'problem6',...
           'Position', pos,...
           'NumberTitle', 'off');
    hold on;
    grid on;
    view(125,30); % CHANGE
    
    title('y[n]  0 <= n <= 50', 'FontSize', 16);
    xlabel('input n', 'FontSize', 16);
    ylabel('imag  output h[n]', 'FontSize', 16);
    zlabel('real  output h[n]', 'FontSize', 16);
    xlim([rng_n(1)-buf_n, rng_n(2)+buf_n]);
    ylim([rng_h_im(1)-buf_im, rng_h_im(2)+buf_im]);
    zlim([rng_h_rl(1)-buf_rl, rng_h_rl(2)+buf_rl]);
    plot3(ax_n     , 0*ax_n   , 0*ax_n   , 'k', 'linewidth', 1);
    plot3(0*ax_h_im, ax_h_im  , 0*ax_h_im  , 'k', 'linewidth', 1);
    plot3(0*ax_h_rl, 0*ax_h_rl, ax_h_rl, 'k', 'linewidth', 1);
    plot3(rng_n(2)+buf_n  , 0                 , 0                 ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0               , rng_h_im(2)+buf_im, 0                 ,'y.', 'markersize', 20, 'linewidth', 10);
    plot3(0               , 0                 , rng_h_rl(2)+buf_rl,'y.', 'markersize', 20, 'linewidth', 10);
    text(rng_n(2)+buf_n, 0                 , 0                 , 'n in', 'FontSize', 16);
    text(0             , rng_h_im(2)+buf_im, 0                 , 'imag out', 'FontSize', 16);
    text(0             , 0                 , rng_h_rl(2)+buf_rl, 'real out', 'FontSize', 16);
    
    stmH = stem3(n_in, h_im, h_rl);
    stmH.LineStyle = '-';
    set(stmH, 'Marker', '^', 'MarkerSize', 5);
    stmH.LineWidth = 2;
    stmH.MarkerFaceColor = 'b';
    stmH.MarkerEdgeColor = 'b';
    stmH.Color = 'r';
    temH = stmH.BaseLine;
    temH.Visible = 'off';  
    %}
end





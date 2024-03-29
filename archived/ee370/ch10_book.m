%{
    chapter 10 : stability, nyquist

    1  :  ex10.1, basic bode
    2  :  ex10.3, 3rd order bode
    3  :  ex10.8, phase and gain margin
    4  :  sa10.5, gm + pm
    5  :  ex10.10, use the bode
    6  :  sa10.6, use the bode for margin
    7  :  sa10.8, nichols

    rlocus()
    rlocusplot()
    sgrid()
    step()
    Control System Designer
    nyquist()  ...right click and say all margins
    nicholas()

    margin()
    APPS --> control --> linear system analyzer
%}
format compact;
clear;
close all;
clc;
select = 7;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms s;
    Gs = 1 / (s + 2);
    tf_sys = f_Fs_2_tf(Gs);
    bode(tf_sys);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms s;
    Gs = (s + 3) / ((s + 2) * (s^2 + 2*s + 25));
    tf_sys = f_Fs_2_tf(Gs);
    pretty(expand(Gs));
    bode(tf_sys);
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms s;
    K = 6;
    Ls = K / ((s + 2) * (s^2 + 2*s + 2));
    KLs = K * Ls;
    tf_KLs = f_Fs_2_tf(KLs);

    syms w;
    Ljw = subs(Ls, s, 1j*w);
    w_r = solve(imag(Ljw)==0, w)
    w_real = w_r(1);
    fprintf("crosses real axis @ w= %0.3f rad\n", w_real);

    real_part = subs(abs(Ljw), w, w_real);
    fprintf("at crossing, mag is:  %0.3f\n", real_part);
    increase = 1 / real_part;
    fprintf("so the gain can be increased:  %0.3fx  -->  %0.3f\n", increase, increase*K);
    G_m = 20 * log10(increase);
    fprintf("so the Gain margin:  %0.3f dB\n", G_m);

    [Gm, Pm, Wcg, Wcp] = margin(tf_KLs)

    %{
    figure;
    hold on;
    bode(tf_KLs, "g-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;

    figure;
    hold on;
    axis equal;
    n = 100;
    theta = linspace(0, 2*pi, n);
    x = cos(theta);
    y = sin(theta);
    nyquist(tf_KLs, "g-");
    plot(x,y, "k:");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;
    %}
end


%------------------------------------------------------------------------------------------
if (select == 4)
    syms s;
    K = 100;
    Ls = K / ((s + 2) * (s + 4) * (s + 6));
    tf_Ls = f_Fs_2_tf(Ls);

    syms w;
    Ljw = subs(Ls, s, 1j*w);
    Ljw_abs = abs(Ljw)
    %{
    [nn, dd] = numden(Ljw)
    nn_r = real(nn)
    nn_i = imag(nn)
    nn_abs = sqrt(nn_r^2 + nn_i^2)
    nn_angle = atan2(nn_i, nn_r)
    dd_r = real(dd)
    dd_i = imag(dd)
    dd_abs = sqrt(dd_r^2 + dd_i^2)
    dd_angle = atan2(dd_i, dd_r)
    Ljw_abs = nn_abs / dd_abs;
    Ljw_angle = nn_angle - dd_angle
    %}
   
    nyquist(tf_Ls)
end


%------------------------------------------------------------------------------------------
if (select == 5)
    syms s;
    syms w;
    K = 200;
    Ls = K / ((s + 2) * (s + 4) * (s + 5));
    Ljw = subs(Ls, s, 1j*w);
    Ljw_abs = abs(Ljw);
    Ljw_angle = angle(Ljw);
    tf_Ls = f_Fs_2_tf(Ls);
    
    %linearSystemAnalyzer(tf_Ls);
    
    real_cross = subs(Ljw_abs, w, 6.16);
    fprintf("\nreal crossing:  %0.3f\n", real_cross);
    
end


%------------------------------------------------------------------------------------------
if (select == 6)
    syms s;
    syms w;
    K = 10000;
    Ls = 1 / ((s + 5) * (s + 20) * (s + 50));
    KLs = K * Ls;
    tf_Ls = f_Fs_2_tf(Ls);
    tf_KLs = f_Fs_2_tf(KLs);

    %
    figure;
    hold on;
    bode(tf_KLs, "r-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;

    figure;
    hold on;
    axis equal;
    nyquist(tf_KLs, "r-");
    set(findall(gcf,'type','line'),'linewidth',2);
    hold off;
    %}

end


%------------------------------------------------------------------------------------------
if (select == 7)
    syms s;
    K = 8000;
    Ls = K / ((s + 5) * (s + 20) * (s + 50));
    tf_Ls = f_Fs_2_tf(Ls);

    figure;
    hold on;
    %axis equal;
    nichols(tf_Ls);
    ngrid;
    hold off;
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end


%%%%%%%%~~~~~~~~END>
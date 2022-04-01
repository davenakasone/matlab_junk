%{
    chapter 2, DSP fundamentals

    1  :  example of sampling
    2  :  zplane()
    3  :  freqz()
    4  :  fft()
    5  :  hw 2.2b
    6  :  hw 2.2a
    7  :  hw 2.2c
    8  :  hw 2.6a
    9  :  hw 2.6b
    10 :  hw 2.6c
    11 :  hw 2.8  ...simple graph/sample, no function calls
    12 :  hw 2.9   ...DFT as fft()
    13 :  hw 2.16a
    14 :  hw 2.16b ... good zplane()
    15 :  hw 2.16c
    16 :  hw 2.17a  ... good freqz()
    17 :  hw 2.17b  
    18 :  hw 2.17c  
    19 :  hw 2.21a
    20 :  hw 2.21b
    21 :  hw 2.21c
    22 :  hw 2.19
    23 :  

    fvtool(),  fft(), fftshift(), freqz(), zplane(), mean(), std(), rand()
%}
close all;
clc;
select = 12;  %        <-----CHANGE

%------------------------------------------------------------------------------------------
if (select == 1)
    n = [0:31]; % Time index n=0,1, . . . ,31
    omega = 0.25*pi; % Digital frequency
    xn = 2*sin(omega*n); % Sine wave generation
    figure('Position', [20, 20, 700, 700]);
    hold on;
    plot(n, xn, "bo", 'markersize', 10, 'linewidth', 3); % Samples are marked by ‘o’
    plot(n, xn, "b-", 'linewidth', 1);
    xlabel("Time index, n");
    ylabel("Amplitude");
    %axis([0 31 -2.5 2.5]); % Define ranges of plot
    %axis([0, 31, -inf, inf]);
    axis padded;
    %save sine.dat xn -ascii; % Save in ASCII data file
end

%------------------------------------------------------------------------------------------
if (select == 2)
    syms z;
    Hz = 1/ (1 - (1/z) + (9/(10*z)));
    figure('Position', [20, 20, 700, 700]);
    zplane(1, [1, -1, 0.9]);
    axis padded;
end

%------------------------------------------------------------------------------------------
if (select == 3)
    bb = [-1];
    aa = [1, -1, 0.9];
    freqz(bb, aa);
end

%------------------------------------------------------------------------------------------
if (select == 4)
    N=100; f = 1000; fs = 10000; % Define parameter values
    n=[0:N-1]; k=[0:N-1];        % Define time and frequency indices
    omega=2*pi*f/fs;             % Frequency of sinewave
    xn=sin(omega*n);             % Generate sinewave
    Xk=fft(xn,N);                % Perform DFT
    magXk=20*log10(abs(Xk));     % Compute magnitude spectrum
    plot(k, magXk); axis([0, N/2, -inf, inf]); % plot from 0 to pi
    xlabel('Frequency index, k');
    ylabel('Magnitude (dB)');
end

%------------------------------------------------------------------------------------------
if (select == 5)
    aa = [1, 2];
    bb = [1, -3/10, -2/5];
    %sys =  filt(aa, bb);
    syms z;
    syms n;
    H_num = 1-(2/z);
    H_den = 1-(3/10)*(1/z)-(2/5)*(1/(z^2));
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    for ii = 0:1:4
        fprintf("\t\t\t\th(%d)=  %0.3f\n", ii, double(subs(Hn, n, ii)));
    end
end

%------------------------------------------------------------------------------------------
if (select == 6)
    syms z;
    syms n;
    H_num = 1;
    H_den = 1-(3/(4*z));
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    for ii = 0:1:4
        fprintf("\t\t\t\th(%d)=  %0.3f\n", ii, double(subs(Hn, n, ii)));
    end
end

%------------------------------------------------------------------------------------------
if (select == 7)
    syms z;
    syms n;
    H_num = 2-(2/z)+(1/(2*z^2));
    H_den = 1;
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    for ii = 0:1:4
        fprintf("\t\t\t\th(%d)=  %0.3f\n", ii, double(subs(Hn, n, ii)));
    end
end

%------------------------------------------------------------------------------------------
if (select == 8)
    syms z;
    H_num = 1;
    H_den = 1-(3/(4*z));
    Hz = H_num / H_den;
    zplane(1, [1, 3/4]);
end

%------------------------------------------------------------------------------------------
if (select == 9)
    syms z;
    H_num = 1 - (2/z);
    H_den = 1 - (3/(10*z) - (2/(5*z^2)));
    Hz = H_num / H_den;
    zplane([1, -2], [1, -3/10, -2/5]);
end

%------------------------------------------------------------------------------------------
if (select == 10)
    syms z;
    H_num = 2 - (2/z) + (1/(2*z^2));
    H_den = 1;
    Hz = H_num / H_den;
    zplane([2, -2, 1/2], 1);
end

%------------------------------------------------------------------------------------------
if (select == 11)
    syms t;
    f_sig = 2000;
    sig_in = sin(2 * pi * f_sig * t);
    
    f_sample = 10000;
    samples = 100;
    nn = 0:(1/f_sample): (1/f_sample)*(samples-1);
    sig_sampled = double(subs(sig_in, t, nn));
    
    dots = 700; 
    bufX = 2*(1/f_sample);
    rngX = [min(nn, [], 'all')-bufX, max(nn, [], 'all')+bufX];
    bufY = .1;
    rngY = [min(sig_sampled, [], 'all')-bufY, max(sig_sampled, [], 'all')+bufY];
    x_ax = linspace(rngX(1), rngX(2), dots);
    y_ax = linspace(rngY(1), rngY(2), dots);
    rider_t = linspace(rngX(1)+bufX, rngX(2)-bufX, dots);
    rider_sig = subs(sig_in, t, rider_t);
    
    figure('Position', [20, 20, 700, 700]);
    hold on;
    grid on;
    view(2);
    tiStr = "sampled sine wave, N = 100";
    title(tiStr, 'fontsize', 26); 
    xlabel('t / n', 'fontsize', 18);
    ylabel('x(t) / x(n)', 'fontsize', 18);   
    xlim([rngX(1), rngX(2)]);
    ylim([rngY(1), rngY(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1); 
    
    plot(rider_t, rider_sig, 'b-', 'linewidth', 2);
    for ii = 1:1:samples
        plot(nn(ii), sig_sampled(ii), 'r.', 'markersize', 20);
    end
end

%------------------------------------------------------------------------------------------
if (select == 12)
    N = 100;
    f = 2000; %2000
    fs = 10000; %10000
    n = 0:1:N-1;
    k = n;
    w = 2*pi*(f/fs);
    xn = sin(w*n);
    Xk = fft(xn, N);
    Xk_mag = 20*log10(abs(Xk));
    plot(k, Xk_mag, 'r', 'linewidth', 1);
    xlabel('frequency index \omega_k', 'fontsize', 18);
    ylabel('magnitude in dB', 'fontsize', 18);  
end

%------------------------------------------------------------------------------------------
if (select == 13)
    [hz, hp, ht] = zplane(1, [1, -3/4]);
    set(findobj(hz, 'Type', 'line'), 'Color', 'r', 'linewidth', 3, 'markersize', 10);
    set(findobj(hp, 'Type', 'line'), 'Color', 'b', 'linewidth', 3, 'markersize', 10);
    set(findobj(ht, 'Type', 'line'), 'Color', 'k', 'linewidth', 1);
end

%------------------------------------------------------------------------------------------
if (select == 14)
    [hz, hp, ht] = zplane([1, 2], [1, -3/10, -2/5]);
    set(findobj(hz, 'Type', 'line'), 'Color', 'r', 'linewidth', 3, 'markersize', 10);
    set(findobj(hp, 'Type', 'line'), 'Color', 'b', 'linewidth', 3, 'markersize', 10);
    set(findobj(ht, 'Type', 'line'), 'Color', 'k', 'linewidth', 1);
end

%------------------------------------------------------------------------------------------
if (select == 15)
    [hz, hp, ht] = zplane([2, -2, 1/2], 1);
    set(findobj(hz, 'Type', 'line'), 'Color', 'r', 'linewidth', 3, 'markersize', 10);
    set(findobj(hp, 'Type', 'line'), 'Color', 'b', 'linewidth', 3, 'markersize', 10);
    set(findobj(ht, 'Type', 'line'), 'Color', 'k', 'linewidth', 1);
end

%------------------------------------------------------------------------------------------
if (select == 16)
    freqz(1, [1, -3/4]);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 17)
    freqz([1, 2], [1, -3/10, -2/5]);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 18)
    freqz([2, -2, 1/2], 1);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 19)
    syms z;
    syms n;
    H_num = 1;
    H_den = 1-(3/(4*z));
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    N = 128;
    NN = 0:1:N-1;
    Hnn = double(subs(Hn, n, NN));
    
    dots = 700; 
    bufX = 1;
    rngX = [min(NN, [], 'all')-bufX, max(NN, [], 'all')+bufX];
    bufY = .2;
    rngY = [min(Hnn, [], 'all')-bufY, max(Hnn, [], 'all')+bufY];
    x_ax = linspace(rngX(1), rngX(2), dots);
    y_ax = linspace(rngY(1), rngY(2), dots);
    
    figure('Position', [20, 20, 900, 500]);
    hold on;
    grid on;
    view(2);
    tiStr = "sampled impulse, N = 0:127";
    title(tiStr, 'fontsize', 26); 
    xlabel('n', 'fontsize', 18);
    ylabel('h(n)', 'fontsize', 18);   
    xlim([rngX(1), rngX(2)]);
    ylim([rngY(1), rngY(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);
    for ii = 1:1:N
        plot(NN(ii), Hnn(ii), 'r.', 'markersize', 10);
    end
end

%------------------------------------------------------------------------------------------
if (select == 20)
    syms z;
    syms n;
    H_num = 1 - (2/z);
    H_den = 1 - (3/(10*z)) - (2/(5*z^2));
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    N = 128;
    NN = 0:1:N-1;
    Hnn = double(subs(Hn, n, NN));
    
    dots = 700; 
    bufX = 1;
    rngX = [min(NN, [], 'all')-bufX, max(NN, [], 'all')+bufX];
    bufY = .2;
    rngY = [min(Hnn, [], 'all')-bufY, max(Hnn, [], 'all')+bufY];
    x_ax = linspace(rngX(1), rngX(2), dots);
    y_ax = linspace(rngY(1), rngY(2), dots);
    
    figure('Position', [20, 20, 900, 500]);
    hold on;
    grid on;
    view(2);
    tiStr = "sampled impulse, N = 0:127";
    title(tiStr, 'fontsize', 26); 
    xlabel('n', 'fontsize', 18);
    ylabel('h(n)', 'fontsize', 18);   
    xlim([rngX(1), rngX(2)]);
    ylim([rngY(1), rngY(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);
    for ii = 1:1:N
        plot(NN(ii), Hnn(ii), 'r.', 'markersize', 10);
    end
end

%------------------------------------------------------------------------------------------
if (select == 21)
    syms z;
    syms n;
    H_num = 1 - (2/z) - (1/(2*z^2));
    H_den = 1;
    Hz = H_num/H_den;
    Hn = iztrans(Hz, z, n);
    pretty(Hn);
    N = 128;
    NN = 0:1:N-1;
    Hnn = double(subs(Hn, n, NN));
    
    dots = 700; 
    bufX = 1;
    rngX = [min(NN, [], 'all')-bufX, max(NN, [], 'all')+bufX];
    bufY = .2;
    rngY = [min(Hnn, [], 'all')-bufY, max(Hnn, [], 'all')+bufY];
    x_ax = linspace(rngX(1), rngX(2), dots);
    y_ax = linspace(rngY(1), rngY(2), dots);
    
    figure('Position', [20, 20, 900, 500]);
    hold on;
    grid on;
    view(2);
    tiStr = "sampled impulse, N = 0:127";
    title(tiStr, 'fontsize', 26); 
    xlabel('n', 'fontsize', 18);
    ylabel('h(n)', 'fontsize', 18);   
    xlim([rngX(1), rngX(2)]);
    ylim([rngY(1), rngY(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);
    for ii = 1:1:N
        plot(NN(ii), Hnn(ii), 'r.', 'markersize', 10);
    end
end

%------------------------------------------------------------------------------------------
if (select == 22)
    f_signal = 1000;
    f_sample = 8000;
    N = 1024;
    nn = 0:1:N-1;
    w_n = 2 * sqrt(0.6) * (rand(1, N) - (1/2));
        check_avg = mean(w_n);    % mean should be about 0
        fprintf("mean of noise:  %f\n", check_avg);
        check_var = (std(w_n))^2; % var should be about 0.2 as given
        fprintf("var of noise:  %f\n", check_var);
    s_n = sin(2 * pi * (f_signal/f_sample) * nn);
    x_n = s_n + w_n;
    X_k = fft(x_n, N);
    Xk_mag = 20*log10(abs(X_k));
    power_w = 0;
    power_sig = 0;
    for ii = 1:1:N
        power_w = power_w + w_n(ii)^2;
        power_sig = power_sig + s_n(ii)^2;
    end
    SNR = 10*log10(power_sig/power_w);
    fprintf("\tSNR:  %f  dB\n", SNR);
    
    dots = 700; 
    bufX = 1;
    rngX = [min(nn, [], 'all')-bufX, max(nn, [], 'all')+bufX];
    bufY = .2;
    rngY = [min(Xk_mag, [], 'all')-bufY, max(Xk_mag, [], 'all')+bufY];
    x_ax = linspace(rngX(1), rngX(2), dots);
    y_ax = linspace(rngY(1), rngY(2), dots);
    
    figure('Position', [20, 20, 900, 500]);
    hold on;
    grid on;
    view(2);
    tiStr = "frequency spectrum, magnitude of fft(), N = 1024";
    title(tiStr, 'fontsize', 26); 
    xlabel('k', 'fontsize', 18);
    ylabel('|X(k)| dB', 'fontsize', 18);   
    xlim([rngX(1), rngX(2)]);
    ylim([rngY(1), rngY(2)]);
    plot(x_ax  , 0*x_ax, 'k', 'linewidth', 1);
    plot(0*y_ax, y_ax  , 'k', 'linewidth', 1);
    
    plot(nn, Xk_mag, 'b-', 'linewidth', 1);
end

%------------------------------------------------------------------------------------------
if (select == 23)
    
    
end
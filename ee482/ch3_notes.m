%{
    chapter 3, FIR filter
    
    1  : ex 3.1, matlab, freqz(), grpdelay()
    2  : ex 3.2, matlab, comb filter ..a multi-band 
    3  : ex 3.7, matlab, low pass filter
    4  : ex 3.8, matlab, gibbs oscillations
    5  : ex hamming, matlab, how to use a hamming window
    6  : ex 3.9, matlab, low pass with hamming
    7  : ex 3.10, matlab, remez() designs FIR filter
    8  : ex 3.12, matlab, dfilt manipulate object for FIR filter
    9  : ex 3.13, matlab, filter and white noise
    10 : ex 3.14, matlab
    11 : ex 3.18, matlab, uses sound file
    12 : ex 3.19, matlab, uses sound file
    13 : ex 3.20, matlab, up filter
    14 : ex 3.21, uses speech file, turn speakers on
    15 : hw2, p3.1   use e#1
    16 : hw2, p3.5  use e#2
    17 : hw2, p3.9
    18 : hw2, p3.11 use #3
    19 : hw2, p3.12
    20 : hw2, p3.13
    21 : hw2, p3.21, see ex3.21, using "TIMIT_4.ASC"

    grpdelay(), freqz(), hamming(), remez(), filter(), filtfilt(), gcd(), interp()
    upfirdn(), decimate()

    wintool, wvtool, fdatool, fvtool
    fir1(), fir2(), kaiserord()
    fir1s(), firpm(), firpmord()
    fircls(), fircls1()
    cfirpm()
    dfilt.dffir
%}
close all;
clc;
select = 21;  %        <-----CHANGE

%------------------------------------------------------------------------------------------
if (select == 1)
    b=[0.5, 0.5]; a = [1];  % Define filter coefficients
    freqz(b,a);             % Show magnitude and phase response
    figure;
    grpdelay(b,a);          % Show group delay
end

%------------------------------------------------------------------------------------------
if (select == 2)
    b=[1 0 0 0 0 0 0 0 -1]; a=[1]; % Define parameters
    freqz(b, a);                   % Show magnitude and phase response
end

%------------------------------------------------------------------------------------------
if (select == 3)
    omegac = 0.4*pi; % Cutoff frequency
    L = 61;		     % Filter order L = 61
    M = (L-1)/2;     % M    since L is odd, L = 2M + 1, solve for M
    l = 0:2*M;       % Coefficient index l
    h = omegac/pi*sinc(omegac*(l-M)/pi);  % Compute coefficients
    omega = -pi:2*pi/200:pi;     % Frequency range
    Hd = freqz(h,1,omega);       % Frequency response 
    plot((omega/pi),abs(Hd)),... % Use normalized frequency
        xlabel('Normalized frequency'), ylabel('Magnitude'), grid;
    axis([-1 1 0 1.2]);
end

%------------------------------------------------------------------------------------------
if (select == 4)
    M = 8;                     % M = 8 (First window)
    L = 2*M+1;                 % Window length
    wn = [ones(1,L)];          % Fectangular window
    omega = -pi:2*pi/200:pi;   % Frequency range
    Hd = freqz(wn,1,omega);    % Frequency response
    mag = 20*log10(abs(Hd));   % Log scale

    M1 = 20;                   % M = 20 (Second window)
    L1 = 2*M1+1;               % Window length
    wn1 = [ones(1,L1)];        % Rectangular window
    omega = -pi:2*pi/200:pi;   % Frequency range
    Hd1 = freqz(wn1,1,omega);  % Frequency response
    mag1 = 20*log10(abs(Hd1)); % Log scale

    subplot(2,1,1), plot((omega/pi),mag), ...
        xlabel('Normalized frequency'), ylabel('Magnitude(dB)'), grid;
    axis([-1 1 -40 40]);
    subplot(2,1,2), plot((omega/pi),mag1), ...
        xlabel('Normalized frequency'), ylabel('Magnitude(dB)'), grid;
    axis([-1 1 -40 40]);
end

%------------------------------------------------------------------------------------------
if (select == 5)
    L = 41; n = 0:40;             % Window length
    wn = hamming(L);              % Generate window coefficients
    omega = -pi:2*pi/200:pi;      % Ffrequency range
    Hd = freqz(wn,1,omega);       % Frequency response
    mag = 20*log10(abs(Hd));      % in dB scale
    subplot(2,1,1), plot(n,wn), ...
    ylabel('Amplitude'); axis([0 (L-1) 0 1]);
    subplot(2,1,2), plot((omega/pi),mag), ...
        xlabel('Normalized frequency'), ylabel('Magnitude(dB)'), grid;
    axis([-1 1 -80 40]);
end

%------------------------------------------------------------------------------------------
if (select == 6)
    omegac = 0.4*pi; % Cutoff frequency
    L = 61;          % Filter order L = 61
    M = (L-1)/2;     % M
    l = 0:2*M;       % Coefficient index l
    h = omegac/pi*sinc(omegac*(l-M)/pi);  % Compute coefficients
    wn = 0.54 -0.46*cos(2*pi*l/(L-1));    % Hamming window
    hwn = h.*wn;     % Windowing
    omega = -pi:2*pi/200:pi;   % Frequency range
    Hr = freqz(h,1,omega);     % Frequency response-Rectangular
    Hd = freqz(hwn,1,omega);   % Frequency response-Hamming  
    plot((omega/pi),abs(Hr),':r',(omega/pi),abs(Hd),'-b');
    xlabel('Normalized frequency'), ylabel('Magnitude'), grid;
    legend('Rectangular window','Hamming window')
    axis([-1 1 0 1.2]);
end

%------------------------------------------------------------------------------------------
if (select == 7)
    f = [0  0.3  0.4  0.6  0.7  1];  % Frequency range
    m = [0  0  1  1  0  0];          % Desired magnitude response
    b = remez(17, f, m);             % Remez FIR filter design
    [h, omega] = freqz(b, 1, 512);   % Frequency response
    plot(f, m, omega/pi, abs(h));
    xlabel('Normalized frequency'); ylabel('Magnitude'), grid;
end

%------------------------------------------------------------------------------------------
if (select == 8)
    b = firls(80,[0 0.11 0.19 1],[1 1 0 0],[1 100]); % Design an FIR filter
    hd = dfilt.dffir(b); % Create the direct-form FIR filter.
    set(hd,'Arithmetic','fixed');  % Quantize filter using 16-bit
    % fvtool(b,hd);      % Compare the fixed-point filter with reference
    h1 = copy(hd);       % Copy hd to h1
    set(h1,'CoeffWordLength',12); % Use 12 bits for coefficients
    fvtool(hd, h1);      % Compare 12-bit & 16-bit filters
end

%------------------------------------------------------------------------------------------
if (select == 9)
    Fs = 8000;      % Sampling rate
    Ts = 1/Fs;      % Total samples
    F = 1500;       % Sinewave frequency
    sinewave = sin(2*pi*F*(0:Ts:1));              % Generate sinewave
    noise = sqrt(0.1).*randn(1,length(0:Ts:1));   % var = 0.1
    xn = sinewave+noise;             % Combine sinewave with noise
    in = xn(1:256);                  % Using first 256 samples only
    xn_int = round(32768*in./max(abs(in))); % Normalize to 16-bit integer
    plot(xn_int); axis([1 256 -32768 32767]);
    ylabel('Amplitude'); xlabel('Sample index, n');
    fid = fopen('xn_int.dat','w');   % Save signal to xn_int.dat
    fprintf(fid,'%4.0f\n',xn_int);   % Save in integer format
    fclose(fid);
end

%------------------------------------------------------------------------------------------
if (select == 10)
    rand('state',0); % Initializing the random number generator
    q = quantizer([16,15],'RoundMode','round'); 
    xq = randquant(q,256,1); % 256 samples in the range [-1,1)
    xin = fi(xq,true,16,15);

    b = firls(80,[0 0.11 0.19 1],[1 1 0 0],[1 100]); % Design FIR filter
    hd = dfilt.dffir(b); % Create the direct-form FIR filter.
    set(hd,'Arithmetic','fixed');  % Quantize filter using 16-bit

    y = filter(hd,xin);  % Fixed-point filtering
    plot(y);
end

%------------------------------------------------------------------------------------------
if (select == 11)
    load wn8kHz.dat -ascii;               % Load the data file
    fftL = 4*16384;
    binsmh=128;

    f=8000*(0:(fftL/2-1))/fftL;           % Frequency scale for display
    y=fft(wn8kHz,fftL);                   % FFT of the specified block
    py=20*log10(abs(y)/fftL);             % Magnitude
    for i = 1:(fftL/2);                   % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));
    grid;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    axis([0 24000 -10 40]);

    x = zeros(size(1:fftL));
    x(1:6:fftL)= wn8kHz(1024:1:(1023+10923));% Decimate to 8 kHz
    h=figure;
    f=48000*(0:(fftL/2-1))/fftL;             % Frequency scale for display
    y=fft(x,fftL);                           % FFT of the specified block
    py=20*log10(2*abs(y)/fftL);              % Magnitude
    for i = 1:(fftL/2);                      % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));
    grid;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    axis([0 24000 -10 40]);

    % Lowpass filtering
    b = fir1(480,1/6);	                  % Design FIR filter in b
    xx = conv(x,b);                       % Decimate to 8 kHz
    h=figure;
    f=48000*(0:(fftL/2-1))/fftL;          % Frequency scale for display
    y=fft(xx,fftL);                       % FFT of the specified block
    py=20*log10(2*abs(y)/fftL);           % Magnitude
    for i = 1:(fftL/2);                   % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));
    grid;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    axis([0 24000 -10 40]);
end

%------------------------------------------------------------------------------------------
if (select == 12)
    load wn48kHz.dat -ascii;           % Load the data from file
    fftL = 16384*4;
    binsmh = 128;

    f=48000*(0:(fftL/2-1))/fftL;       % Frequency scale for display
    y=fft(wn48kHz,fftL);               % FFT of the specified block
    py=20*log10(abs(y)/fftL);          % Magnitude spectrum
    for i = 1:(fftL/2);                % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));grid;title('(a) Original singal spectrum');
    xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
    axis([0 24000 10 40]);

    fftL = fftL/4; 
    len= size(wn48kHz);
    x = wn48kHz(1:6:len);              % Decimate to 8 kHz
    h=figure;
    f=8000*(0:(fftL/2-1))/fftL;        % Frequency scale for display
    y=fft(x,fftL);                     % FFT of the specified block
    py=20*log10(3*abs(y)/(2*fftL));    % Magnitude spectrum
    for i = 1:(fftL/2);                % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));
    grid;title('(b) Decimation by 6 without lowpass filter');
    xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
    axis([0 4000 10 40]);

    xx = decimate(wn48kHz,6,100,'fir'); % Decimate to 8 kHz
    h=figure;
    f=8000*(0:(fftL/2-1))/fftL;         % Frequency scale for display
    y=fft(xx,fftL);                     % FFT of the specified block
    py=20*log10(3*abs(y)/(2*fftL));     % Magnitude spectrum
    for i = 1:(fftL/2);                 % Smooth diplay
    psy(i) = sum(py(i:(i+binsmh)))/binsmh;
    end
    plot(f,psy(1:(fftL/2)));
    grid;title('(c) Decimation by 6 with lowpass filter');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    axis([0 4000 10 40]);
end

%------------------------------------------------------------------------------------------
if (select == 13)
    g = gcd(48000, 44100);              % Greatest common divisor, g = 300
    U = 44100/g;                        % Upsampling factor, U=147
    D = 48000/g;                        % Downsampling factor, D = 160

    N = 24*D;
    b = fir1(N,1/D,kaiser(N+1,7.8562)); % Design FIR filter in b
    b = U*b;                            % Passband gain = U
    Fs = 48000;                         % Original sampling frequency: 48kHz
    n = 0:10239;                        % 10240 samples, 0.213 seconds long
    x = sin(2*pi*1000/Fs*n);            % Original signal, sinusoid at 1 kHz
    y = upfirdn(x,b,U,D);               % 9408 samples, still 0.213 seconds

    % Overlay original (48 kHz) with resampled signal (44.1 kHz) in red
    stem(n(1:49)/Fs,x(1:49)); 
    hold on
    stem(n(1:45)/(Fs*U/D),y(13:57),'r','filled');
    xlabel('Time (seconds)');
    ylabel('Signal value');
end

%------------------------------------------------------------------------------------------
if (select == 14)
    load timit_4.asc -ascii;    % Load speech file
    soundsc(timit_4, 16000)     % Play at 16 kHz
    disp('Press a key to continue ...');
    pause;
    timit2 = decimate(timit_4,8,60,'fir'); % decimate by 8
    soundsc(timit2, 2000)       % Play the decimated speech
    disp('Press a key to continue ...');
    pause;
    timit48 = interp(timit_4,3);% Interpolate to 48 kHz
    soundsc(timit48,48000);     % Play the interpolate speech
    disp('Press a key to continue ...');
    pause;
    disp('Example completed');
end

%------------------------------------------------------------------------------------------
if (select == 15)
    syms w;
    f_sample = 8000;
    Hw_mag = sqrt((1/2) * (1 + cos(w)));
    Hw_mag_db = 20 * log10(Hw_mag);
    omg = solve(Hw_mag_db == -3, w);
    F = omg(2) / 2;
    fprintf("\n\t|H(w)| = -3 dB @ F = +/-  %0.3f\n", F);
    fprintf("\n\tw = %f rad/samp\n", omg(2));
    
    b=[0.5, 0.5]; a = 1;  
    freqz(b,a);           
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 16)
    L = 8;
    H_num = [1, 0, 0, 0, 0, 0, 0, 0, -1];
    H_den = 1;
    roots(H_num)
    [hz, hp, ht] = zplane(H_num, H_den);
    set(findobj(hz, 'Type', 'line'), 'Color', 'r', 'linewidth', 3, 'markersize', 10);
    set(findobj(hp, 'Type', 'line'), 'Color', 'b', 'linewidth', 3, 'markersize', 10);
    set(findobj(ht, 'Type', 'line'), 'Color', 'k', 'linewidth', 1);
    figure;
    freqz(H_num, H_den); 
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 17)
    syms z;
    N = 11;
    hn_a = [-4, 1, -1, -2, 5, 0, -5, 2, 1, -1, 4];
    hn_b = [-4, 1, -1, -2, 5, 6, 5, -2, -1, 1, - 4];
    figure;
    freqz(hn_a, 1); 
    my_lines = findall(gcf, 'type', 'line')
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
    figure;
    freqz(hn_b, 1); 
    my_lines = findall(gcf, 'type', 'line')
    set(my_lines(1), 'Color', 'r', 'linewidth', 2, 'LineStyle', '--');
    set(my_lines(2), 'Color', 'b', 'linewidth', 2, 'LineStyle', '--');
end

%------------------------------------------------------------------------------------------
if (select == 18)
    M_lo = 32;
    M_hi = 64;
    l_lo = 0 : 1 : (2 * M_lo - 1);
    l_hi = 0 : 1 : (2 * M_hi - 1);
    wc = 0.6 * pi;
    h_lo = fir1(M_lo, wc/pi, 'high', ones(1, M_lo + 1)); % no window effect
    h_hi = fir1(M_hi, wc/pi, 'high', ones(1, M_hi + 1)); % no window effect
    
    freqz(h_lo, 1);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
    figure;
    freqz(h_hi, 1); 
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2, 'LineStyle', ':');
    set(my_lines(2), 'Color', 'b', 'linewidth', 2, 'LineStyle', ':');
end

%------------------------------------------------------------------------------------------
if (select == 19)
    M_lo = 32;
    M_hi = 64;
    wc = 0.6 * pi;
    %h_lo = fir1(M_lo, wc/pi, 'high'); % hamming, by default
    %h_hi = fir1(M_hi, wc/pi, 'high'); % hamming, by default
    h_lo = fir1(M_lo, wc/pi, 'high', blackman(M_lo+1)); % blackman
    h_hi = fir1(M_hi, wc/pi, 'high', blackman(M_hi+1)); % blackman
    
    freqz(h_lo, 1);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
    figure;
    freqz(h_hi, 1); 
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2, 'LineStyle', ':');
    set(my_lines(2), 'Color', 'b', 'linewidth', 2, 'LineStyle', ':');
end

%------------------------------------------------------------------------------------------
if (select == 20)
    M = 20;
    w_left = 0.4;
    w_right = 0.5;
    h = fir1(M, [w_left, w_right]);
    
    freqz(h, 1);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
end

%------------------------------------------------------------------------------------------
if (select == 21)
    M = 100;
    L = 2 * M;
    wc_up = pi / L;
    wc_down = pi / M;
    load timit_4.asc -ascii;                 % load sound file
    timit60 = interp(timit_4, 3);            % up-sample to 60 kHz
    %f1 = fir1(M, timit60, wc_up, 'low');
    soundsc(timit60, 60000)                  % play at 60 kHz
    disp('Press a key to continue ...');
    pause;
   
    timit12 = decimate(timit60, 5);    % down-sample 60 kHz/5= 12 kHz
    soundsc(timit12, 12000);           % play at 12 kHz
    disp('Press a key to continue ...');
    pause;
    
    timit5 = decimate(timit60, 12);    % down-sample 60 kHz/12= 5 kHz
    soundsc(timit5, 5000);           % play at 5 kHz
    disp('Press a key to continue ...');
    pause;
    
    timit3 = decimate(timit60, 20);    % down-sample 60 kHz/20= 3 kHz
    soundsc(timit3, 3000);           % play at 3 kHz
    disp('Press a key to continue ...');
    pause;
    
    disp('Example completed');
end

%------------------------------------------------------------------------------------------
if (select == 0)
    load timit_4.asc -ascii;    % Load speech file
    soundsc(timit_4, 16000)     % Play at 16 kHz
    disp('Press a key to continue ...');
    pause;
    timit2 = decimate(timit_4,8,60,'fir'); % decimate by 8
    soundsc(timit2, 2000)       % Play the decimated speech
    disp('Press a key to continue ...');
    pause;
    timit48 = interp(timit_4,3);% Interpolate to 48 kHz
    soundsc(timit48,48000);     % Play the interpolate speech
    disp('Press a key to continue ...');
    pause;
    disp('Example completed');
    
    
    
    
    no_window_lo = zeros(1, M_lo + 1);
    no_window_hi = zeros(1, M_hi + 1);
    h_lo = fir1(M_lo, wc/pi, 'high');
    h_hi = fir1(M_hi, wc/pi, 'high');
    
    freqz(h_lo, 1);
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2);
    set(my_lines(2), 'Color', 'b', 'linewidth', 2);
    figure;
    freqz(h_hi, 1); 
    my_lines = findall(gcf, 'type', 'line');
    set(my_lines(1), 'Color', 'r', 'linewidth', 2, 'LineStyle', ':');
    set(my_lines(2), 'Color', 'b', 'linewidth', 2, 'LineStyle', ':');
end

%%%%>>>>END
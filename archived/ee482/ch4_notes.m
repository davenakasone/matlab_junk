%{
    chapter 4, IIR filter
    
    1  : ex 4.3, matlab, butter(), bilinear(),
    2  : ex 4.7, matlab, tf2zp(),
    3  : ex 4.8, matlab, buttord(), butter(),
    4  : ex 4.9, matlab, same as above, 
    5  : ex 4.15, matlab, butter and filtering
    6  : ex 4.16, matlab, 
    7  : ex 4.17, matlab, fvtool(), 
    8  : ex 4.18, matlab, fvtool(),
    9  : hw, p4.5, given IIR filter tf, find zeros/poles, see if stable
    10 : hw, p4.7
    11 : hw, p4.9
    12 : hw, p4.10
    13 : hw, p4.12, 13, 14
    14 : hw, p4.18
    15: play
    
    bilinear()
    butter()
    buttord()
    cfirpm()
    decimate()
    dfilt.dffir
    ellip()
    fdatool
    fft()
    fftshift()
    filter()
    filterDesigner
    filtfilt()
    fir1()
    fir1s()
    fir2()
    fircls()
    fircls1()
    firpm()
    firpmord()
    freqz()
    fvtool
    gcd()
    grpdelay()
    hamming()
    iirpeak()
    interp()
    kaiserord()
    mean()
    rand()
    remez()
    residuez()
    sosfilt()
    std()
    tf2zp()
    upfirdn() 
    wintool
    wvtool
    zplane()
    
%}
close all;
clc;
select = 15;  %        <-----CHANGE


%------------------------------------------------------------------------------------------
if (select == 1)
    Fs = 2000;                          % Sampling frequency
    Wn = 300;                           % Edge frequency
    n = 4;                              % Order of analog filter
    [b, a] = butter(n, Wn, 's');        % Design analog filter
    [bz, az] = bilinear(b, a, Fs);      % Determine digital filter
    freqz(bz,az,512,Fs);                % Display magnitude & phase
end


%------------------------------------------------------------------------------------------
if (select == 2)
    b = [0.5, 0, -0.18];     % Numerator
    a = [1, 0.1, -0.72];     % Denominator
    [z, p, c] = tf2zp(b, a)  % Display zeros, poles, and gain
end


%------------------------------------------------------------------------------------------
if (select == 3)
    Wp = 800/4000;  % Passband frequency
    Ws= 1600/4000;  % Stopband frequency
    Rp = 1.0;       % Passband ripple
    Rs = 20.0;      % Stopband attenuation
    
    [N, Wn] = buttord(Wp, Ws, Rp, Rs); % First step
    [b, a] = butter(N, Wn);            % Second step
    freqz(b, a, 512, 8000);            % Display frequency response
end


%------------------------------------------------------------------------------------------
if (select == 4)
    Wp = [100  200]/500;    % Passband frequencies
    Ws = [50  250]/500;     % Stopband frequencies
    Rp = 3;                 % Passband ripple
    Rs = 30;                % Stopband attenuation
    [N, Wn] = buttord(Wp, Ws, Rp, Rs); % Estimate filter order
    [b, a] = butter(N, Wn);	% Design a Butterworth filter
    fvtool(b, a);			% Analyze the designed IIR filter
end


%------------------------------------------------------------------------------------------
if (select == 5)
    Fs = 1000;                     % Samplng rate
    f0 = 150;
    N=300; 
    A=sqrt(2); 
    w0=2*pi*f0/Fs;                 
    n = [0:N-1];                   % Time index
    sn = A*sin(w0*n);              % Sine sequence
    sd = 12357; rng(sd);           % Define seed value
    vn = (rand(1,N)-0.5)*sqrt(12); % Zero-mean, unit-variance white noise
    xn = sn+vn;                    % Sinewave embedded in white noise

    %  Design a bandpass filter
    Wp = [140  160]/(Fs/2);        % Passband edge frequencies
    Ws = [120  180]/(Fs/2);        % Stopband edge frequencies
    Rp = 3;                        % Passband ripple
    Rs = 40;                       % Stopband ripple
    [N, Wn] = buttord(Wp, Ws, Rp, Rs); % Find the order
    [b, a] = butter(N, Wn);        % Design an IIR filter

    y = filter(b, a, xn);          % IIR filtering
    subplot(2,1,1); 
    plot(n, xn);
    xlabel('Time index, n'); ylabel('Amplitude');
    subplot(2,1,2); 
    plot(n, y);
    xlabel('Time index, n'); ylabel('Amplitude');
end


%------------------------------------------------------------------------------------------
if (select == 6)
    Fs = 1000; 
    f0 = 150;
    N=300;                        
    A=sqrt(2); 
    w0=2*pi*f0/Fs; 
    n = [0:N-1];                   % Time index
    sn = A*sin(w0*n);              % Sine sequence
    vn = (rand(1,N)-0.5)*sqrt(12); % Zero-mean, unit-variance white noise
    xn = sn+vn;                    % Sinewave embedded in white noise

    %  Design a bandpass filter
    Wp = [140  160]/(Fs/2);        % Passband edge frequencies
    Ws = [130  170]/(Fs/2);        % Stopband edge frequencies

    Rp = 3;  
    Rs = 40;
    [N, Wn] = buttord(Wp, Ws, Rp, Rs);
    [b, a] = butter(N, Wn);

    %  Convert it to sos
    sos = tf2sos(b,a);
    y = sosfilt(sos,xn);           % IIR filtering
    subplot(2,1,1); 
    plot(n, xn);
    xlabel('Time index, n');
    ylabel('Amplitude');
    subplot(2,1,2); 
    plot(n, y);
    xlabel('Time index, n');
    ylabel('Amplitude');
end


%------------------------------------------------------------------------------------------
if (select == 7)
    Fs = 10000;                % Sampling rate
    Wo = 1000/(Fs/2);          % First filter peak frequency
    BW = 500/(Fs/2);           % First filter bandwidth
    W1 = 2500/(Fs/2);          % Second filter peak frequency
    BW1 = 200/(Fs/2);          % Second filter bandwidth
    [b,a] = iirpeak(Wo,BW);    % Design first filter
    [b1,a1] = iirpeak(W1,BW1); % Design second filter
    fvtool(b,a,b1,a1);         % Analyze both filter
end


%------------------------------------------------------------------------------------------
if (select == 8)
    Fs = 10000;             % Sampling rate
    f0 = 1500;              % Frequency in Hz
    w0 = 2*pi*f0/Fs;
    rz=0.8; rp=0.9;         % Define parameters
    b=[1, -2*rz*cos(w0), rz*rz]; % Define numerator coefficints
    a=[1, -2*rp*cos(w0), rp*rp]; % Define denominator coefficients
    fvtool(b,a);            % Analyze the filter
end


%------------------------------------------------------------------------------------------
if (select == 9)
    syms z;
    H_num_a = 1 + 1.414*z^-1 + z^-2;
    H_num_b = 1 + 2*z^-1 + z^-2;
    H_num = H_num_a * H_num_b;
    H_den_a = 1 - 0.8*z^-1 + 0.64*z^-2;
    H_den_b = 1 - 1.0833*z^-1 + 0.25*z^-2;
    H_den = H_den_a * H_den_b;
    Hz = H_num / H_den;
    %pretty(expand(H_num));
    pretty(expand(H_den));
    %pretty(simplify(Hz));
    Hz_rootz = double(solve(H_den == 0, z))
    H_den_a_rtz = solve(H_den_a==0, z);
    p_a1 = H_den_a_rtz(1,1);
    p_a2 = H_den_a_rtz(2,1);
    check1a = double(abs(p_a1*p_a2))
end


%------------------------------------------------------------------------------------------
if (select == 10)
    syms a;
    Hz = ( -a + (1/z) ) / ( 1 - a*(1/z) );
    figure;
    hold on;
    %freqz([-1, 1], [1, -1]); % [-a, 1] [1, -a]
    %freqz([-5, 1], [1, -5]); % [-a, 1] [1, -a]
    freqz([-100, 1], [1, -100]); % [-a, 1] [1, -a]
end


%------------------------------------------------------------------------------------------
if (select == 11)
    syms z;
    Hz_n = 12 - 2/z + 3/z^2 + 20/z^4;
    Hz_d = 6 - 12/z + 11/z^2 -5/z^3 + 1/z^4;
    Hz_check = Hz_n / Hz_d;
    
    H_num = [12, -2, 3, 0, 20];
    %H_num_r = roots(H_num)
    H_den = [6, -12, 11, -5, 1];
    %H_den_r = roots(H_den)
    [rr, pp, kk] = residuez(H_num, H_den)
    %
    [hz, hp, ht] = zplane(H_num, H_den);
    set(findobj(hz, 'Type', 'line'), 'Color', 'r', 'linewidth', 3, 'markersize', 10);
    set(findobj(hp, 'Type', 'line'), 'Color', 'b', 'linewidth', 3, 'markersize', 10);
    set(findobj(ht, 'Type', 'line'), 'Color', 'k', 'linewidth', 1);
    %
    Hz_split = kk + ( rr(1,1)/(1-pp(1,1)/z) ) + ( rr(2,1)/(1-pp(2,1)/z) ) +...
        ( rr(3,1)/(1-pp(3,1)/z) ) + ( rr(4,1)/(1-pp(4,1)/z) );
    check = Hz_split - Hz_check;
    checkk = double(subs(check, z, 1)) % about 0
    pretty(expand(Hz_split))
end


%------------------------------------------------------------------------------------------
if (select == 12)
    f_s = 8000;
    f_pass = 1600;
    f_stop = 2000;
    pb_rip = 0.5;
    sb_atten = 40;
    nn = 6;
    
    [bb_s, aa_s] = ellip(nn, pb_rip, sb_atten, 2*f_pass/f_s)
    [bb_z, aa_z] = bilinear(bb_s, aa_s, f_s, f_stop);
    freqz(bb_s, aa_s, 512, f_s);
    fvtool(bb_s, aa_s);
end


%------------------------------------------------------------------------------------------
if (select == 13)
    %fdatool
    filterDesigner
end


%------------------------------------------------------------------------------------------
if (select == 14)
    Hz_num = 0.0662 .* [1, 3, 3, 1];
    Hz_den = [1, -0.9356, 0.5671, -0.1016];
    %fvtool(Hz_num, Hz_den);
    
    syms z;
    syms n;
    Hzn = 0.0662 * ( 1 + 3/z + 3/z^2 + 1/z^3 );
    Hzd = 1 - 0.9356/z + 0.5671/z^2 - 0.1016/z^3;
    Hz = Hzn / Hzd;
    pretty(Hz);
    h_n = iztrans(Hz, z, n);
    pretty(simplify(h_n));
    nn = 0:1:99;
    h_nn = double(subs(h_n, n, nn));
    hold on;
    grid on;
    axis padded;
    plot(nn, h_nn, "r-", "linewidth", 3);
end


%------------------------------------------------------------------------------------------
if (select == 15)
    syms z;
    Hz = (4-6/z+2/z^2) / (1 + 1/(4*z) -1/(8*z^2));
    pretty(Hz);
    pretty(simplify(expand(Hz)))
    D = (1 + 1/(4*z) -1/(8*z^2));
    pretty(simplify(D))
    solve(D==0, z)
    N = (4-6/z+2/z^2);
    solve(N==0, z)
    
    [rr, pp, kk] = residuez([4, -6, 1], [1, 1/4, -1/8])
end
%%%%~~~~END>

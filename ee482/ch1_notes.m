%{
    chapter 1 basics, DSP intro
    
    1  : check your 2D graph, CT
    2  : sample the CT function
    3  : quantizing % TODO

%}
close all;
clc;
select = 3;


%------------------------------------------------------------------------------------------
if (select == 1)
    syms t;
    f_1 = sin(1 * t);
    f_2 = sin(5 * t);
    start = 0;
    stop = 3*pi;
    g_ct_2d([f_1, f_2], t, [start, stop; start-pi, stop+pi]);
end


%------------------------------------------------------------------------------------------
if (select == 2)
    syms t;
    sig_1 = sin(1 * t);
    sig_2 = sin(5 * t);
    start = 0;
    stop = 3*pi;
    f_sample_1 = 1; % Hz
    f_sample_2 = 1; % Hz
    g_dt_sampler([sig_1, sig_2], t, [start, stop; start, stop], [f_sample_1, f_sample_2], 1);
end


%------------------------------------------------------------------------------------------
if (select == 3)
    syms t;
    myFun1 = 5 - exp(.5*t);
    myFun2 = 4 - exp(.5*t);
    start = 0;
    stop = 6;
    
    g_dt_samplerQ(funz, t, rangz, freqz, Qbits, show_trace)
end


%END
clear;
close all;
clc;
% spectrum test
data = readtable('plot_header.csv');
var_ = data.Var2;

var_fft = fft(var_)

figure(Position=[20, 20, 800, 800]);
spectrogram(var_fft)

%{
N = 1024;
n = 0:N-1;
w0 = 2*pi/5;
x = sin(w0*n)+10*sin(2*w0*n);
s = spectrogram(x);
spectrogram(x,'yaxis')
%}

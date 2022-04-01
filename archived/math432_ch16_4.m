%{
    ch16.4  applications of residue integration

    #1      ex4   simple pole on real axis
%}
clc;
close all;
clearvars;
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
old_val = sympref('HeavisideAtOrigin', 1);
format short; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]

global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
%global Zxy; syms Zxy; Zxy = X + 1j*Y;% just use it as a hold and imply  when needed
global Z; syms Z;
global alpha; syms alpha; assume(alpha, 'real');
global beta; syms beta; assume(beta, 'real');

                            select = 1;  % CHANGE CHANGE CHANGE


%------------------------------------------------------------------------------------------ #1
if select == 1
    
end
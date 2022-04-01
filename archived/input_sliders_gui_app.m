%{
    use the app desinger if you start having a really big project
    this will work if you just have a few functions and want to see a few changes


    sys is a good way to make a transfer function     H  = N / D
%}

clc;
close;
clearvars;
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice

%{
prompt = ' (1) rectangular , (2) cylindrical , (3) spherical  : ';
response = input(prompt, 's');
if response == '1'                 % takes as a string
    display([response,' : retangular is stupid']);
end
if response == '2'                 % takes as a string
    display([response,' : cylindrical is stupid']);
end
if response == '3'                 % takes as a string
    display([response,' : spherical is is stupid']);
end
pause(3);
clc;
pause(3);
response = input(prompt);
if response == 1                 % takes as a string
    fprintf("\n %d : rectangular is dumb\n", response);
end
if response == 2                 % takes as a string
    fprintf("\n %d : cylindrical is dumb\n", response);
end
if response == 3                 % takes as a string
    fprintf("\n %d : spherical is dumb\n", response);
end
%}

zeta = .5;  % damping ratio
wn = 2;     % natural frequency
sys = tf(wn^2, [1, 2*zeta*wn, wn^2]);  % establish transfer function

%create GUI to show step response  and configure axes
f = figure;
ax = axes('Parent', f, 'position', [.13, .39, .77, .54]);
h = stepplot(ax, sys);
setoptions(h, 'XLim',[0,10],'YLim',[0,2]);
%add slider
b = uicontrol('Parent', f, 'Style', 'slider', 'Position',[81, 54, 419, 23], 'value', zeta, 'min', 0, 'max', 1);
bgcolor = f.Color;
bl1 = uicontrol('Parent', f, 'Style', 'text', 'Position', [50, 54, 23, 23], 'String', '0', 'BackgroundColor', bgcolor);
bl2 = uicontrol('Parent', f, 'Style', 'text', 'Position', [500, 54, 23, 23], 'String', '1', 'BackgroundColor', bgcolor);
bl3 = uicontrol('Parent', f, 'Style', 'text', 'Position', [240, 25, 100, 23], 'String', 'Damping Ratio', 'Backgroundcolor', bgcolor);
b.Callback = @(es, ed)updateSystem(h, tf(wn^2,[1,2*(es.Value)*wn,wn^2]));
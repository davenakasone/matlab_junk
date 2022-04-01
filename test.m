%{
    home > preferences : shortcuts, line wraps, font size...

    Specify OpenGL Implementation for Current Session
    To specify an OpenGL implementation for the current session of MATLAB, use one of these techniques.

    Software OpenGL — Start MATLAB from the command prompt on your system using the command matlab -softwareopengl. 
    This command works only on Windows and Linux systems. Macintosh systems do not support software OpenGL.

    Basic hardware-accelerated OpenGL — Type opengl hardwarebasic at the MATLAB command prompt.

    Hardware-accelerated OpenGL — Type opengl hardware at the MATLAB command prompt.

    Specify OpenGL Implementation for Future Sessions
    To set your preferences so that MATLAB always starts with the specified implementation of OpenGL, use one of these techniques.

    Software OpenGL — Type opengl('save','software') at the MATLAB command prompt. Then, restart MATLAB.

    Basic hardware-accelerated OpenGL — Type opengl('save','hardwarebasic') at the MATLAB command prompt. 
    Then, restart MATLAB.

    Hardware-accelerated OpenGL — Type opengl('save','hardware') at the MATLAB command prompt. Then, restart MATLAB.

    Undo preference setting — Execute opengl('save','none') at the MATLAB command line. Then, restart MATLAB.
%}
close all;
clear all;
clc;
select = 1;
%opengl('save','software');
%opengl('save','hardwarebasic')
%info = rendererinfo(gca)
%opengl('save','hardware')

%-----------------------------------------------------------------------------------------------------
if select == 1
    fprintf("\n\tit works\n");
    xx = [0, 1, 2, 3, 4, 5];
    yy = [2, 3, 4, 5, 6, 7];
    figure('Position',[20, 20, 800, 800]);
    plot(xx, yy, 'g-', 'LineWidth', 3);
end


%-----------------------------------------------------------------------------------------------------
if select == 99

end
%%%%%%%%~~~~~~~~END>  xxx.m
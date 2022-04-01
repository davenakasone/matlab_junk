function result = fun_harminic2test(funIn)
%{
    result = 0 if not harmonic, 1 if harmonic
    funIn is seperated f(z) = u(x,y) + j v(x,y)
    test1
        namb2 (u ) = u_xx + u_yy = 0
    test2
        namb2 (v) = v_xx + v_yy = 0

    if it passes RC test, it should pass this   
        ... real  и imag    parts of analytic funcitons are harmonic
%}
result = 0;
global X; % need X and Y as the symbols to z = x + jy
global Y;
global Z;

    u = real(funIn);
    u_x = diff(u, X, 1);
    u_xx = diff(u_x, X, 1);
    u_y = diff(u, Y, 1);
    u_yy = diff(u_y, Y, 1);
    
    v = imag(funIn);
    v_x = diff(v, X, 1);
    v_xx = diff(v_x, X, 1);
    v_y = diff(v, Y, 1);
    v_yy = diff(v_y, Y, 1);
    
    test1 = u_xx + u_yy;
    test2 = v_xx + v_yy;
    
    if test1 == 0 && test2 == 0
        fprintf('f(z) = %s   IS HARMONIC\n', funIn);
        fprintf('u_xx + u_yy = %s   +    %s  = %.1f\n', u_xx, u_yy, test1);
        fprintf('v_xx + v_yy = %s   +    %s  = %.1f\n', v_xx, v_yy, test2);
        result = 1;
    else
        fprintf('f(z) = %s   NOT harmonic\n', funIn);
        fprintf('u_xx + u_yy = %s   +    %s  = %.1f\n', u_xx, u_yy, test1);
        fprintf('v_xx + v_yy = %s   +    %s  = %.1f\n', v_xx, v_yy, test2);
        fprintf('the namblas do not = 0\n');
    end
end


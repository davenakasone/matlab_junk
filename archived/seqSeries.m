%{
    most useful seq will be one that converges to a limit

    #1      1 +  1/(1+j)^n     good use of numbers for points
    #2      p_n = j^(p_n-1)  n>= 2, p_1 = j      j^j = exp(j log(j)) = exp( pi /2)
    #3      exp(z) = lim, n-->oo, (1 + z/n)^n   ... best way, get seq first, then use max
    #4      exp(z) by taylor     cumsum()
    #5      1/cosh(z)    as a power series  and TAYLOR

%}
clc;
clf;
close all;
clearvars;
                            sel = 5;

%------------------------------------------------------------------------------------------ #1
if sel == 1
    N = 15;
    p = ones(1,N);
    rp = real(p);    % not necessary, but nice to get separation...always good to establish array
    ip = imag(p);
    
    for n = 1:N
        p(n) = 1 + 1/(1 + 1j)^n;
        rp(n) = real(p(n));
        ip(n) = imag(p(n));
        b = num2str(n);
        plot(rp(n), ip(n));
            text(rp(n), ip(n), b, 'FontSize', 6);
            hold on;
            pause(.25);
            n
            p(n)
    end
    grid
        plot(p,'-');
        title('first 15 terms 1 + 1/(1+j)^n');
        hold off;
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    N = 50;
    p = ones(1,N);
    rp = real(p);
    ip = imag(p);
    p(1) = 1j;
    
    for n = 2:N
        p(n) = 1j^p(n-1);
        rp(n) = real(p(n));
        ip(n) = imag(p(n));
        b = num2str(n);
        plot(rp(n), ip(n), '*');
        pause(.3);
        hold on;
        text(rp(n), ip(n), b, 'FontSize', 6);
        hold on;
        n
        p(n)
    end
    
    plot(p(1), '*');
    grid on;
    title('seq p_n = j^(p_n-1)');
    hold off;
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    maxP = 3;
    N = 1000;
    seq = zeros(1,N);
    z = 1 + 1j;
    
    figure();
    hold on;
    grid on;
    xlim([0, maxP]);
    ylim([0, maxP]);
    xlabel('real');
    ylabel('imaginary');
    title('e as a seq');
    rax = linspace(0, maxP , N);
    iax = linspace(0 , maxP , N);
    plot(rax  , 0*rax, 'k', 'linewidth', 1);
    plot(0*iax, iax  ,'k', 'linewidth', 1);
    
    plot(real(exp(z)), imag(exp(z)), 'b*');
    for n = 1:N
        seq(1,n) = ( 1 + z/n)^n;
         b = num2str(n);
        plot( real(seq(1,n)), imag(seq(1,n)), 'r.');                  % one or
        %pause(.1);
        hold on;
        pnt = num2str(n);
        %text( real(seq(1,n)), imag(seq(1,n)), pnt, 'FontSize', 6);    % the other
        hold on;
    end
    hold off;
    d = min(real(seq))
    dd = max(real(seq))
    e = min(imag(seq))
    ee = max(imag(seq))
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    Nmax = 11;
    a = 0:Nmax-1;
    q = 1 ./ factorial(a);
    b = (2 + 2j).^a;
    SN = q .* b;
    format long
    S = (cumsum(SN)).'; % sums all elements in row vector       . stops conj()
    N = [1:Nmax]';           % turns into coulumb vector
    fprintf('terms of series:\n');
    [N, S]
    fprintf(' exp( 2 + j2) : %s\n', exp(2+2j) );
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    % for 1/cosh(z)  there will be holes at z = +/- j pi/2    so abs(z) < pi/2 to convg
    % center it at z0 = 0  "mclaurin" of the taylor
    % 1/cosh(z) = SUM c_n z^n  [0,oo]
    %{
    clf;
    x = linspace(-2,2,100);
    z = x + 1j*0;
    w = 1 ./ cosh(z);
    u = real(w);
    v = imag(w);
    plot(x,u);
    grid on;
    hold on;
    ser = 1 - (z.^2 / 2) + ( 5 * z.^4 / 4 );  % 3 term approx
    plot(x, ser, '--r', 'LineWidth' , 2);
    ylim([-5, 10]);
    xlabel('x');
    ylabel('1/cosh(z)  with approximation');
    %}
    % this would be a pain in the ass to get each taylor term
    % if you wanted the 4th coeff with z0 = 0
    %{
    syms a;
    syms z;
    f = 1/cosh(z);
    a = diff(f, z, 4);
    %or
    syms b;
    b(z) = @(z) ( diff(f, z, 4) / factorial(4) )
    % or   b = @(z) ...
    coeff4 = b(0);  % 5/24 = 1/4!  d4{z} 1/cosh(z)
    %}
    % or do it the real way and use 'Taylor'  just symbols
    %{
    syms z;
    syms T;
    f = 1/cosh(z);
    center = 1j;
    T = taylor(f, z, center,'order', 3);  % radius pi/2 - 1   ceneter j
    %}
    % TAYLOR for numberics
    syms z;
    syms b;
    fun = 1/cosh(z);
    center = 1j;
    b(z) =@(z) taylor(fun, z, center, 'order', 3);
    b(.1 + .1j);
    eval( b(.1 + .1j) );
end
    
    
        
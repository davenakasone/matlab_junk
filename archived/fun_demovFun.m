function rootz = fun_demovFun()
%{
    MANUAL Funciton...  must change n, r, th

    this one is for functions use fun_demovNum() for numbers

    funIn formated  Z^( 1 / Nroot ) - const       
    
        implies    Z^( 1 / Nroot ) = const

        rootsOut returns [ Nroot x 2 ] solution of the roots
%}

    % ie      z^4 - 1.1 = 0   --->   z^4 = 1.1        -->   z = 1.1^(1/4)
    
    z = -8; %CHANGE
    n = 3;                   %CHANGE     % n = 4; 
    r = abs(z);          % r = abs(1.1 + 1j*0);
    th = angle(z);       % th = angle(1.1 + 1j*0);
    rootz = zeros(n,1);             % no change
    
    
    for idx = 1:n
        k = idx-1;
        realPart = r^(1/n) *  cos( (th/n) + (2*sym(pi)*k)/n );
        imagPart = 1j*r^(1/n) * sin((th/n) + (2*sym(pi)*k)/n );
        rootz(idx) = realPart + imagPart;  % k = 0, 1, 2, ...., n-1
        %fprintf('check %s  -->  %s\n', rootz(idx), rootz(idx)^n);
    end
end


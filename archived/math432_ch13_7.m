%{
        ch 13.7 complex natural logarithm

        #1      ex1   numbers and infinite solutions  
        #2      graphs   ...branch cut
        #3      p5   
        #4      p6
        #5      p7
        #6      p8
        #7      p9
        #8      p10
        #9      p11
        #10     p12
        #11     p13
        #14     p14
        #15     p15
        #16     p16
        #17     p17
        #18     p18
        #19     p19
        #20     p20
        #21     p21
        #22     p22
                
%}
clc;
close all;
clearvars;


                sel = 22;  % CHANGE CHANGE CHANGE
                
                
global X; syms X; assume(X, 'real'); 
global Y; syms Y; assume(Y, 'real');
global Z; Z = (X + 1j*Y);

%------------------------------------------------------------------------------------------ #1
if sel == 1
    x = 1;
    y = 0;
    z = x + 1j*y;
    fprintf('ln ( 1 ) = %.3f +/- 2n pi j   n = 1,2,...    prinipal Ln(1) = 0 \n', log(z));
    
    x = 4;
    y = 0;
    z = x + 1j*y;
    fprintf('\nln(4) =  %.3f  +/-  +/- 2n pi j   n = 1,2,...   prinipal is first...\n',log(z));
    
    x = -4;
    y = 0;
    z = x + 1j*y;
    fprintf('\nln(-4) = %s    +/- (2n+1) pi j\n', log(z));
    
    x = 4;
    y = -3;
    z = x + 1j*y;
    ans = log(z);
    fprintf('\nlog(4-j3) :  %.3f   -j%.3f\n', real(ans), abs(imag(ans)));  % best practice
end


%------------------------------------------------------------------------------------------ #2
if sel == 2
    %fun_graphSurCon(log(Z), 5);
    %fun_tetraView(log(Z), 0, 0)
    %fun_graphSurCon(log(Z), 1); % you can see regular ln(x)
    %fun_graphSurCon(log(Z), 2); % branch cut 
    % several branch cuts:
    fun_graphSurCon2(log(Z),log(Z) + 2*sym(pi), 2);
end


%------------------------------------------------------------------------------------------ #3
if sel == 3
    % ln(-11)
    x = -11;
    y = 0;
    z = x + 1j*y;            % just plug it in
    temp1 = log( abs(z) );
    temp2 = 1j * angle(z);
    sol = temp1 + temp2
    log(z)  
end


%------------------------------------------------------------------------------------------ #4
if sel == 4
    % ln( 4 +  j4)
    x = 4;
    y = 4;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z) 
end


%------------------------------------------------------------------------------------------ #5
if sel == 5
    % ln( 4 -  j4)
    x = 4;
    y = -4;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z) 
end


%------------------------------------------------------------------------------------------ #6
if sel == 6
    % ln( 1 +  j)
    x = 1;
    y = 1;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
    % ln( 1 -  j)
    x = 1;
    y = -1;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
end


%------------------------------------------------------------------------------------------ #7
if sel == 7
    % ln( .6 +  j .8)
    x = .6;
    y = .8;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
end


%------------------------------------------------------------------------------------------ #8
if sel == 8
    % ln( -15 +  j .1)
    x = -15;
    y = .1;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
    % ln( -15 -  j .1)
    x = -15;
    y = -.1;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
end


%------------------------------------------------------------------------------------------ #9
if sel == 9
    % ln( exp(1) * j )
    x = 0;
    y = exp(1);
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2
    log(z)
end

%------------------------------------------------------------------------------------------ #10
if sel == 10
    % ln(exp(1))
    x = exp(1);
    y = 0;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2; % principal
    pts = zeros(11);
    pts(6) = sol;
    for k = 1:5
        pts(k) = sol + 2*pi*1j*k;
    end
    for k = 7:11
        pts(k) = sol - 2*pi*1j*(k-6);
    end
    fun_graph2Dcomplex(pts);
    % check
    fprintf('should get back %.3f  j %.3f\n\n', real(z), imag(z));
    for k = 1:11
        fprintf(' %.3f   j %.3f\n',  real(exp(pts(k))) , imag(exp(pts(k)))   );
    end
end


%------------------------------------------------------------------------------------------ #11
if sel == 11
    % ln(1)
    x = 1;
    y = 0;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2; % principal
    pts = zeros(11);
    pts(1) = sol;
    for k = 2:6
        pts(k) = sol + 2*pi*1j*k;
    end
    for k = 7:11
        pts(k) = sol - 2*pi*1j*(k-6);
    end
    fun_graph2Dcomplex(pts);
    % check
    fprintf('should get back %.3f  j %.3f\n\n', real(z), imag(z));
    for k = 1:11
        fprintf(' %.3f   j %.3f\n',  real(exp(pts(k))) , imag(exp(pts(k)))   );
    end
end


%------------------------------------------------------------------------------------------ #14
if sel == 14
    % ln(-7)
    x = -7;
    y = 0;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2; % principal
    pts = zeros(11);
    pts(6) = sol;
    for k = 1:5
        pts(k) = sol + 2*pi*1j*k;
    end
    for k = 7:11
        pts(k) = sol - 2*pi*1j*(k-6);
    end
    fun_graph2Dcomplex(pts);
    % check
    fprintf('should get back %.3f  j %.3f\n\n', real(z), imag(z));
    for k = 1:11
        fprintf(' %.3f   j %.3f\n',  real(exp(pts(k))) , imag(exp(pts(k)))   );
    end
end


%------------------------------------------------------------------------------------------ #15
if sel == 15
    % ln( exp(j) )
    x = cos(1);
    y = sin(1);
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2; % principal
    pts = zeros(11);
    pts(1) = sol;
    for k = 2:6
        pts(k) = sol + 2*pi*1j*(k-1);
    end
    for k = 7:11
        pts(k) = sol - 2*pi*1j*(k-6);
    end
    fun_graph2Dcomplex(pts);
    % check
    fprintf('should get back %.3f  j %.3f\n\n', real(z), imag(z));
    for k = 1:11
        fprintf(' %.3f   j %.3f\n',  real(exp(pts(k))) , imag(exp(pts(k)))   );
    end
end


%------------------------------------------------------------------------------------------ #16
if sel == 16
    % ln( 4 + j3 )
    x = 4;
    y = 3;
    z = x + 1j*y;            % just plug it in
    temp1 = sym(log( sym(abs(z)) ));
    temp2 = sym(1j * angle(z));
    sol = temp1 + temp2; % principal
    pts = zeros(11);
    pts(1) = sol;
    for k = 2:6
        pts(k) = sol + 2*pi*1j*(k-1);
    end
    for k = 7:11
        pts(k) = sol - 2*pi*1j*(k-6);
    end
    fun_graph2Dcomplex(pts);
    % check
    fprintf('should get back %.3f  j %.3f\n\n', real(z), imag(z));
    for k = 1:11
        fprintf(' %.3f   j %.3f\n',  real(exp(pts(k))) , imag(exp(pts(k)))   );
    end
end

%------------------------------------------------------------------------------------------ #17
if sel == 17
    % ln( j^2 )  ~= 2 ln(j)    holy shit
    z1 = -1
    z2 = 2*log(1j)
end


%------------------------------------------------------------------------------------------ #18
if sel == 18
    % solve log(z) = -pi j / 2
    lhs = -sym(pi)*1j/2;
    fprintf('log(z) = %s\n',lhs);
    temp = exp(lhs);
    fprintf(' z = %s\n', temp);
end


%------------------------------------------------------------------------------------------ #19
if sel == 19
    % solve log(z) = 4-j3
    lhs = 4-1j*3
    temp = exp(lhs)
end
    


%------------------------------------------------------------------------------------------ #20
if sel == 20
    % ln(z) = e - j pi
    rhs = exp(1) - 1j*pi
    temp = exp(rhs)
    answ=double(rhs - log(temp))
end


%------------------------------------------------------------------------------------------ #21
if sel == 21
    % ln(z) = .6 + .4j
    rhs = .6 + 1j*.4
    temp = exp(rhs)
    answ=double(rhs - log(temp))
end


%------------------------------------------------------------------------------------------ #22
if sel == 22
    % solving  2j^2j     just use z^c = exp ( c ln z )
    z = 2*1j;
    c = 2*1j;
    ln = log(abs(z)) + 1j*angle(z); % natural log part
    pow = c * ln;
    sol = exp(pow)
    answ = z^c
end
    
    
    
function result = fun_rCRtest(funIn)   % rectangular inputs only   z = X + jY  or f(z) in

result = 1; % 1 if passes test, 0 if fails test
global X; % need X and Y as the symbols to z = x + jy
global Y;
%global Z;
%Z = X + 1j*Y;
%syms k;
u = real(funIn); 
v = imag(funIn);

dudx = simplify(diff(u, X, 1));
dvdy = simplify(diff(v, Y, 1));
fprintf('dudx = %s  ,  dvdy = %s \n', dudx, dvdy);
test1 = dudx - dvdy;
if test1 ~= 0
    result = 0; 
    fprintf('failed first test, dudx ~= dvdy\n');
end

dudy = simplify(diff(u, Y, 1));
dvdx = simplify(diff(v, X, 1));
test2 = dudy + dvdx;
fprintf('dudy = %s  ,  dvdx = %s \n', dudy, dvdx);
if test2 ~= 0
    result = 0; 
    fprintf('failed second test, dudy ~= -dvdx\n');
end

if result == 1
    fprintf('f(z) = f(x + jy) = %s    PASSED Cauchy-Riemann test\n', funIn);
else
    fprintf('f(z) = f(x + jy) = %s    FAILED Cauchy-Riemann test\n', funIn);
end

end


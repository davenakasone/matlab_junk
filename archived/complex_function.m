%{
    a complex function and some roots
%}


syms z;  % z is a symbol ... should maintian distinction
f(z) = 2*z^2 + 2*z + 2;  % CHANGE

sols = solve(f==0,z);  % changes as roots increase...see degree
display(sols);


num_elms = numel(sols); % number of roots
modulus = abs(double(sols(1)));

for i = 2:num_elms
    temp = abs(double(sols(i)));
    if modulus <= temp
        modulus = temp;
    end
end

display(num_elms);
display(modulus);

x_u = -modulus-2:.1:modulus+2;
y_u = -modulus-2:.1:modulus+2;

[X_u, Y_u] = meshgrid(x_u,y_u);  % a fancy command for complex plane

f_z = 2*((X_u +1j*Y_u).^2) + 2*(X_u +1j*Y_u)+2; % redefined in terms of complex components    "." for matrix
surfc(x_u, y_u, abs(f_z));
xlabel('real');
ylabel('imaginary');
colormap jet; % sexy

num = 1 + 1j;
anum = (num)^(37);
display(anum);

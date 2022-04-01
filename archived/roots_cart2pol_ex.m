% cartesian/rectangular to polar   cart2pol() and pol2cart()
% atan() is bad

z = 1 + 1*j;
[th, r] = cart2pol(real(z),imag(z));
display(r);
display(th);

ang_z = angle(z);
mag_z = abs(z);
display(mag_z);
display(ang_z);

z_new = pol2cart(ang_z, mag_z);
%z_new = mag_z*cis*ang_z;
z_new = mag_z*(cos(ang_z)+1*j*sin(ang_z)); % best to use Euler
display(z_new);

% for demovire, use roots()
%  w^6 - (1 + j) = 0
coeffs = [1, 0, 0, 0 ,0, 0, -(1+j*1)];
answ = roots(coeffs)
display(answ);
display(answ.^6); % to check

% no guaruntee on oops in matlab    must check all
% for (3+4j)^(5/7)  take it by root, then raise it
cuffs = zeros(1,8); %prep from an array/row vector with 0's to feed to root function
m=5;
n=7; % should be 7 roots, exponent gives it away
cuffs(1) = 1; % means w^7 term has coeff of 1
cuffs(8) = -(3+j*4); % setting up w^7 = (3 +4j)
poly_roots = roots(cuffs); % should give you 7 roots
final = poly_roots.^5; % raise to 5 to complete z^(m/n)
display(final);
check = final.^(7/5);
display(check); % only one will actually recover what you want, need to see all




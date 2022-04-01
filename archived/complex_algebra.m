% syms is more for variables
% sym is more for consts like pi    stops round off error

selector = 4;

if selector == 1
syms x real;  % x can only be real
syms y;  % y can be real or imaginary

syms a b c d e f real;
w = (1/(a+j*b)) + ((c+j*d)/(e+j*f));     % w, u, v are not symbols, only equal to symbols
u = real(w);
display(u);
display(expand(u));
pretty(u); % best

% input numbers
u = subs(u, [a e], [3 7]);  % put in a = 3, b = 7  ...subs is how you get values to symbols
pretty(u);
end

if selector == 2
    a = sqrt(5);
    display(a);
    % vs
    b = sym(sqrt(5));
    display(b);
    c = b^7;
    display(c);
    display(eval(c));  % eval() is how you get the number....don't double()
end

if selector == 3
    z1 = ((5)^(1/2)+j*(7)^(1/2))^(11); % ewww
    display(z1);
    a = sym(sqrt(5));
    b = sym(sqrt(7));
    z2 = (a + j*b)^(11); % sym cleans it up
    display(z2);
    display(eval(z2)); % to check
    display(who); % use who to see all the shit in your work space...or look at the pane
end

if selector == 4
    display(sin(1000*pi)); % false
    display(sin(1000*sym(pi))); % takes actual value, no rounding error
end
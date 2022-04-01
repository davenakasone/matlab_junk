syms x;
syms y(x);
ode = diff(y,x)  == (2*y-x)/(4*x-12*y);

cond = y(7/2) == -2/3;

sol(x)=dsolve(ode,cond);
display(sol);
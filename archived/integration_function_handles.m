%{
    function handle to pass infomration to another function    or:

    Specifying callback functions (for example, a callback that responds 
        to a UI event or interacts with data acquisition hardware).

    Constructing handles to functions defined inline instead of stored in a 
        program file (anonymous functions).

    int() function integrates symbolic expressions   
    integral() function is for numeric functions, usually anonymous functions. 

%}

clc;

%{ 
    Named function handles represent functions in existing program files, 
    including functions that are part of MATLAB and functions that you create 
    using the function keyword. To create a handle to a named function, precede 
    the function name with @.
%}
f = @sin;
m = fminbnd(f,0,2*pi);
fprintf('the value "m" that minimize sin(x) [0,2pi] is : %.1f°\n', rad2deg(m));
fprintf(' check, sin(m) = sin( %.1f°) = %.1f\n\n', rad2deg(m), sin(m));
dump = input('any key:', 's');
clc;

%{
    Anonymous function handles (often called anonymous functions) represent single inline 
    executable expressions that return one output. To define an anonymous function, 
    enclose input argument names in parentheses immediately after the @ operator, 
    and then specify the executable expression.
        Anonymous functions can accept multiple inputs but return only one output
%}
f = @(x,y) sqrt(x.^2 + y.^2);  % helps to use unit operations
radius = f(3,0);
fprintf(' the radius is: %.1f\n', radius);
ptx = 4;
pty = 3;
if f(ptx, pty) ~= 3
    fprintf('point is not on circle with radius = 3, radius here is: %.1f\n',...
        f(ptx, pty) );
end
dump = input('any key:', 's');
clc;

%{
    intead of making a seperate funciton in another file like:
    function y = cubic(x)
        y = x.^3 ; 
    end

    then doing   q = integral(@cubic,0,1)    to integrate [0,1]

    use this:
%}

fun = @(x) (x.^3);
%area = int(fun,0,1);  no symbol, int won't work like this
area = integral(fun,0,1) % works
% or try
syms t;
area = int(fun(t), t, 0, 1)




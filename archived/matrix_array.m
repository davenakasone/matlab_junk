%{
    matlab is matrix/array based
    the syntax is blasphemy

    lots of built in operators and functions .' magic() \ ect

    keep up to date:
    https://www.mathworks.com/help/matlab/referencelist.html?type=function&category=index&s_tid=CRUX_lftnav_function_index

    indexing has end, lots of : operator, and doesn't have to be contiguous...starts at 1 not 0
    assign [] to delete, but try to delete entire column or row

    linear index if you want...stored contigous in memory

    logical index is helpful

    [a, a] horizontal concat
    [a; a] vertical concat
%}

a = [1, 2, 3, 4];
n = size(a);       %returns dimensions
display(a);
display(n(2));

b = [1, 2, 3; 4, 5, 6; 7, 8, 9];
display(b);
display(size(b));

c = 1:10;      % c = 1 to 10
display(c);
display(c(3)); % index off position

d = 1:2:10;  % d = 1 to 10 in steps of 2
display(d);

e = 5:-1:1;  % e = 5 to 1 in steps of -1
display(e);

% linspace() is like operator : , but you control points...matlab guesses rest, should take begin, end, and even space
f = linspace(1, 20, 7) % 7 points in 1-20
display(f);

%transpose operator '
f = f';
display(f);


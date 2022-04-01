%{
    no allocations
        not friendly with dynamic vars
    hard to break   no goto / jump

%}

clc; 
clf;
close all;
clearvars;
select = 2; % CHANGE HERE
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice


mtx = ones(1,10);          % everything is a matrix...2D spec is good
temp = zeros(1,10);
count = 1;
display(mtx);
display(temp);

if select == 1
    while count <= 10;
        temp(1, count) = mtx(1, count);    
        %temp(count) = mtx(count);% can use linear indexing and just spec depth
        count = count + 1;
    end
    display(temp);
end


if select == 2               % while loop does less crazy shit
    for i = 1:3:10      % start:increment:stop     or just use start:stop
        temp(i) = mtx(i);
    end
    display(temp);
end
    
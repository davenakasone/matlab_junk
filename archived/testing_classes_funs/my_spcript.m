%{

    don't really have to import in matlab

    fun_feed()  takes an input and returns the input with 1 added to it
    fun_blank() trying to make a blank function call

    pause() , isa(), class(),  should help out as things get more complicated

    if you have old code, or a hold-out that doesn't like classes,
        take the struct, make a class obj (instantiate), and use obj1 = target
            target will put struct in class

    the class will make an array of targets as you start instatiating them
    classes protect you from adding adition fields to structures and other silly stuff
    no use having struct when you have properties in the class

    varargin is a helpful word  "variable arguments input"
    can even inherit from matlab classes, but that might not be smart
%}

clear;
%clear classes;    do if you need to force remaining definitions out
clc;

x = 2;
y = fun_feed(x);
display('wait 1 second')
pause(1);
clc;
display(y);
fprintf(' var "y" is a double? %d\n\n', isa(y, 'double'));

fun_blank();  % not inputs needed, just a call


mullah_mike = cls_terrorist(2);     % not enogh arguments...default constructor used
mullah_mike.name = 'osama';         % obj.method_name()
mullah_mike.identify1();
display(mullah_mike);

big_slick = cls_terrorist('mullah omar', 'iraq', 3, 2324); % constructor used
display(big_slick)
display(big_slick.bounty_body);     % need property there before dependent can work
big_slick.dispID();   % a protected property
id_range = big_slick.give_hint();
display(id_range);    % it passes through a protected method

jihadi_john = cls_super_terrorist('ak47', 'john', 'canada', 27, 321);  % order is important
display(jihadi_john);
loc = jihadi_john.cords();

perp{1} = cls_terrorist('big$', 'dubai', 34, 99);          % load up an array of objs
perp{2} = cls_terrorist('rick_ross', 'mexico', 77, 777);
perp{3} = cls_super_terrorist('pablo','columbia', 11, 111);

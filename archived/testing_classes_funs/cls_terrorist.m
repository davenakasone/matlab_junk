%{
    this is pathetic compared to C ++
    they think they used the key oop principles    blasphemy at its finest

    properties are data....
        basically a struct

    methods are actions, only for instance of obj of class

    without a default constructor, you will just get an instance with blank values

    the target is for construtor use...creating with default properties
    overload default constructor to take input arguments
            constructor takes same name as class

    getters and setters are good to start going public, private, protected


    you can't add/remove, but you can change all you want
        control access level with attributes

    properties         just a naked struct
    properties (access = protected)                   can see, but can't change
    properties (constnat)                         use for your speed of light and all
    properties (Dependent)
    methods                          free for all
    methods (Access = protected)

    if there is no specification, or public, user can see
    if protected, user can't see
    this is there answer to encapsulation  (hide details, only let outside see what you want)
        ...clean interface to use obj

    classes can be used as a building block
        inherit, then add what you want for new one
            deriving is also possible by composition

            classdef moving_target < target     %  moving_target inherits from target
                    moving_target is a subset of target
                        then just define properties and methods that are needed
%}

classdef cls_terrorist
    
    properties      % all public by default
        name;
        country;
        fbi_rank;
        body_count;
    end
    
    properties (Access = protected)
        jwics_id = 666;                % has a getter to display...
                                            % control through functions as needed
    end
    
    properties (Constant)
        bounty = 50000;    % in dollars       can't be changed
    end
    
    properties (Dependent)        % gets calculated everytime value is asked for
        bounty_body;                  % needs a getter function
    end
    
    methods
        % constructor does not need variable to represent obj...name it as class
        function this_target = cls_terrorist(nm,cnt, rank, bodies)
           % you don't get unlimited constructors like C ++
           % so:
            if nargin == 4     % only do work if you get all 4 arguments
                this_target.name = nm;
                this_target.country = cnt;
                this_target.fbi_rank = rank;
                this_target.body_count = bodies;
            end
        end
        
        % first input is reference to class obj
        % notice function name goes with no output argument if it is doing 1-way op
        function identify1(this_target)
            display([' terroist name: ', this_target.name]);
        end
        
        function bounty_body = get.bounty_body(this_obj)      % getter   has to be public
             bounty_body = cls_terrorist.bounty / this_obj.body_count;
             % the const property is specific to the class ( could use this_obj.light_spd also
        end
        
        function dispID(this_obj)
            display(this_obj.jwics_id);
        end
        
        function answer = give_hint(this)
            answer = this.hint();
        end
    end
    
        
    methods (Access = protected)   % keep to class
        function answer = hint(this)
            if this.jwics_id < 100
                answer = "the jwics id is less than 100";
            else
                answer = "the jwics id is greater than or = 100";
            end
        end
    end
    
end


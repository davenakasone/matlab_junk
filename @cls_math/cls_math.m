classdef cls_math < handle
%{
    !!! matlab treats any function in the class folder as a method of the class !!!
    "printing" argument :  0, don't display steps: 1, display steps
    
    the work flow is just to feed a function in terms of "Z" and calculate
    static methods aren't too much of a focus right now
    
    testing
    #999_1: test_self(this, num_in)                 
    #999_2: result = test_self_ex(this, num_in)
    #999_3: test_static()                    
    #999_4: result = test_static_ex(num_in) 
    
    self
    #0 constructor
    #2 feed(this, fun_new)     update the obj's function use "Z"
    
    static
    #1 fun_image_z(tran, zformat, rng_rx, rng_ty) 
    
%}
%------------------------------------------------------------------------------------------
%******************************************************************************************
%------------------------------------------------------------------------------------------

    properties
        name;     % only used in construction
        funz; % a complex function in terms of z
    end
    
%------------------------------------------------------------------------------------------
%******************************************************************************************
%------------------------------------------------------------------------------------------
    methods

        function obj = cls_math(obj_name) % constructor
            if nargin == 1
                obj.name = obj_name;
            end
        end
%------------------------------------------------------------------------------------------#0

        function result = test_self(this, num_in) % function/method part of classdef
            result = num_in + 999;
            fprintf('%s says %d\n', this.name, result);
        end
%------------------------------------------------------------------------------------------#999_1
        
        result = test_self_ex(this, num_in) % argument to output, function in different file
%------------------------------------------------------------------------------------------#999_2

        function feed(this, fun_new)
            this.funz = 0; % clear old function (if any)
            this.funz = fun_new; % input the new one
        end
%------------------------------------------------------------------------------------------#2
        
        function print_fun(this)
            pretty(this.funz);
        end
%------------------------------------------------------------------------------------------#4
    end
    
    
%------------------------------------------------------------------------------------------
%******************************************************************************************
%------------------------------------------------------------------------------------------


    methods(Static)
        
        function test_static() % no output, no arguments, function part of classdef
            fprintf('static methods are under construction\n');
        end
%------------------------------------------------------------------------------------------#999_3
        
        result = test_static_ex(num_in) % argument to output, function in different file
%------------------------------------------------------------------------------------------#999_4

       fun_image_z(tran, zformat, rng_rx, rng_ty) % map "Z" = x+jy, r exp(1j th) to range
%------------------------------------------------------------------------------------------#1

       lhrh = fun_lft(zs, ws) % [z1,z2,z3], [w1,w2,w3] get lhs, rhs, and w(z) as LFT
%------------------------------------------------------------------------------------------#3

        

    end
%------------------------------------------------------------------------------------------
%******************************************************************************************
%------------------------------------------------------------------------------------------
end




%******************************************************************************************
%******************************************************************************************
%******************************************************************************************
%******************************************************************************************
%******************************************************************************************
%******************************************************************************************

%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------

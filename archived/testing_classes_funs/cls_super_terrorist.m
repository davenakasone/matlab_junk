classdef cls_super_terrorist < cls_terrorist    % inherits
    % all the properties and methods come in
    
    
    properties
        weapon;
    end
    
    properties (Access = protected)
        direction = "kill or capture";
        lat = 33.999;
        long = 36.232;
    end
    
    methods
        % overload the constructor to get to base class
        function newObj = cls_super_terrorist(wep, varargin)
            newObj = newObj@cls_terrorist(varargin{:});            % old class
            if nargin == 1   % use default if no argument provided
                newObj.weapon = wep;  
            end
        end
            
        function location = cords(this)
            location(1,1) = this.lat;
            location(1,2) = this.long;
        end
    end
end


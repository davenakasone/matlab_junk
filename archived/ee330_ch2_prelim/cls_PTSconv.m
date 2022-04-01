classdef cls_PTSconv < handle    % not a value class anymore...fuuck you very much
    % bridge between MATLABS cart2pol, pol2cart, cart2sph, sph2cart
    % MATLAB is good, but:
    %                       clynd is azimuth [ -pi to pi ]    not [0, 2pi] like book
    %                       sphr is azimuth  [ -pi to pi ]    not [0, 2pi] like book
    %                       sphr is elv  [ -pi/2 to pi/2 ]    not [0, pi] like book
    % only for points...no vectors
    
    %{
        fyi, the setter even more of an abortio than the getter
        classdef your_class < handle    for value classes that need to be changed by refernc
    
    
    
            inType     R, C, S    that is it 
    
      pts:
            rec_x     rec_y     rec_z
            cyn_r     cyn_t     cyn_z                 feed it in the book way also 
            sph_r     sph_e     sph_a
    
        constructor
        feed( 1st, 2nd, 3rd, type )
        print()
        clear()
        magD()
    
        don't touch other shit           
    %}
    
    properties
        pts;      % once fed, tell it input format, other 2 coordinates populate
        inType;   % stores string specifiny system input coordinates are provided in
    end
    
    
    methods
        
        function this = cls_PTSconv(str)
            if nargin == 1
                %fprintf('%s', str);
            end
            this.pts = nan(3, 3);
            %fprintf('\npoint converting obj made\n');
            %disp(this.pts);
        end
        
        
        function print(this)
            fprintf('\n    *** angles are actually held in radians *** \n');
            fprintf('rec (x, y, z)   :  %.4f  ,  %.4f   ,  %.4f  \n',...
                this.pts(1, 1), this.pts(1, 2), this.pts(1, 3) );
            fprintf('cyn ( r, fi, z) :  %.4f  ,  %.4f°  ,  %.4f  \n',...
                this.pts(2, 1), rad2deg(this.pts(2, 2)), this.pts(2, 3) );
            fprintf('sph (r, th, fi) :  %.4f  ,  %.4f°  ,  %.4f° \n',...
                this.pts(3, 1), rad2deg(this.pts(3, 2)), rad2deg(this.pts(3, 3)) );
        end
        
        function magD(this)
            fprintf('\n    *** all the magnitudes should be the same *** \n');
            fprintf('\nrec length :  %.4f \n',...
                sqrt( sum(this.pts(1, :).^2) ) ); % the clean way
            fprintf('cyn length :  %.4f \n',...
                sqrt( this.pts(2, 1)^2 + this.pts(2, 2)^2 + this.pts(2, 3)^2 ));
            fprintf('sph length :  %.4f \n',...
                sqrt( this.pts(3, 1)^2 + this.pts(3, 2)^2 + this.pts(3, 3)^2 ));
        end
        
        
        function clear(this)
            global tester;
            %display(subs(tester, 5));
            this.pts = nan(3, 3);
            this.inType = 'a';
        end
        
        
        function feed(this,cord1, cord2, cord3, type)
            if  type == 'R'
                this.inType = 'rec';
                this.pts(1, 1) = cord1;
                this.pts(1, 2) = cord2;
                this.pts(1, 3) = cord3;
                this.rec_to_cyn();
                this.rec_to_sph();
            end
 
            if  type == 'C'
                this.inType = 'cyn';
                this.pts(2, 1) = cord1;
                this.pts(2, 2) = cord2;
                this.pts(2, 3) = cord3;
                this.cyn_to_rec();
                this.rec_to_sph();  % ok to call after others, using rec as a bridge
            end
            
            if  type == 'S'
                this.inType = 'sph';
                this.pts(3, 1) = cord1;
                this.pts(3, 2) = cord2;
                this.pts(3, 3) = cord3;
                this.sph_to_rec();
                this.rec_to_cyn(); % circular call
            end 
        end
 
        
        function rec_to_cyn(this)
            rx = this.pts(1, 1);
            ry = this.pts(1, 2);
            rz = this.pts(1, 3);
            [cf, cr, cz] = cart2pol(rx, ry, rz);
            this.pts(2, 3) = cz;
            this.pts(2, 1) = cr;
            if cf >= 0
                this.pts(2, 2) = cf;
            else
                this.pts(2, 2) = cf + 2*pi;
            end
        end
        
        
        function rec_to_sph(this)
            rx = this.pts(1, 1);
            ry = this.pts(1, 2);
            rz = this.pts(1, 3);
            [sf, st, sr] = cart2sph(rx, ry, rz);
            if sf >= 0
                this.pts(3, 3) = sf;
            else
                this.pts(3, 3) = sf + 2*pi;
            end
            this.pts(3, 2) = (pi/2) - st;
            this.pts(3, 1) = sr;    
        end
        
        
        function cyn_to_rec(this)
            cr = this.pts(2, 1);
            cf = this.pts(2, 2);
            cz = this.pts(2, 3);
            if cf > pi
                cf = cf - 2*pi;
            end
            [rx, ry, rz] = pol2cart(cf, cr, cz);
            this.pts(1, 1) = rx;
            this.pts(1, 2) = ry;
            this.pts(1, 3) = rz;
        end
        
        
        function sph_to_rec(this)
            sr = this.pts(3, 1);
            st = this.pts(3, 2);
            sf = this.pts(3, 3);
            if sf > pi
                sf = sf - 2*pi;
            end
            st = (pi/2) - st;
            [rx, ry, rz] = sph2cart(sf, st, sr);
            this.pts(1, 1) = rx;
            this.pts(1, 2) = ry;
            this.pts(1, 3) = rz; 
        end
        
    end
end


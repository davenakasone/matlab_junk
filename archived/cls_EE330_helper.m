classdef cls_EE330_helper < handle    % not a value class anymore...handle class that can ref self
    %{ 
        user's responsibility to tranpose() , simplify(), double() ect
        bridge between MATLAB's cart2pol, pol2cart, cart2sph, sph2cart
        MATLAB is good, but:
                            cyn is azimuth  [ -pi to pi ]        not   [0, 2pi] like book
                            sph is azimuth  [ -pi to pi ]        not   [0, 2pi] like book
                            sph is elv      [ -pi/2 to pi/2 ]    not   [0, pi] like book
     only for points...no vectors
        
    don't touch constructor, helpers, ect...can't protect because it's a class handle
    wf1: use to convert points, take conversions for inputs
    wf2: use to get a transformation vector, multiply it by an equation in your problem
    wf3; use to make inputs of an equation entirley in terms of one coord sys
    wf4: use for quick proofs
      
    primary function: feed()        
    
      pts:
            rec_x     rec_y     rec_z
            cyn_r     cyn_f     cyn_z               this class is desinged to agree with book             
            sph_r     sph_t     sph_f
    
        obj.feed(point, type);      % point as 1x3 mtx,  type 'R' , 'C' , or 'S' 
        *** don't touch the helpers ***
        obj.print();                % prints 'pts' in 3x3 ...take from pts    recP = obj.pts(1,:);
        obj.clear();                % clears 'pts' to NaN , automatic on every feed
        obj.graph();                % graphs current point (if any), connector, cyn check, sph check
        mag();                      % more validation to make sure point is same distance from origin
        
        trans = transVecRC(this, vecIn)    % transform vec [rx, ry, rz] --> [cr, cf, cz]
        trans = transVecRS(this, vecIn)    % transform vec [rx, ry, rz] --> [sr, st, sf]
        trans = transVecCR(this, vecIn)    % transform vec [cr, cf, cz] --> [rx, ry, rz]
        trans = transVecCS(this, vecIn)    % transform vec [cr, cf, cz] --> [sr, st, sf]
        trans = transVecSR(this, vecIn)    % transform vec [sr, st, sf] --> [rx, ry, rz]
        trans = transVecSC(this, vecIn)    % transform vec [sr, st, sf] --> [cr, cf, cz]
    
                {static methods} to assist in calculations
        
            RC   you have something in terms of [rx, ry, rz] and want to put it in [cr, cf, cz]
        eqn = obj.simpleRCsub();    % provide anything in terms of [rx, ry, rz] , get with [cr, cf, cz]
        eqn = obj.getU_RC();        % provides the unit vector conversion matrix
        
            RS   you have something in terms of [rx, ry, rz] and want to put it in [sr, st, sf]
        eqn = obj.simpleRSsub();    % provide anything in terms of [rx, ry, rz] , get with [sr, st, sf]
        eqn = obj.getU_RS();        % provides the unit vector conversion matrix
        
            CR   you have something in terms of [cr, cf, cz] and want to put it in [rx, ry, rz]
        eqn = obj.simpleCRsub();    % provide anything in terms of [cr, cf, cz] , get with [rx, ry, rz]
        eqn = obj.getU_CR();        % provides the unit vector conversion matrix
        
            CS   you have something in terms of [cr, cf, cz] and want to put it in [sr, st, sf]
        eqn = obj.simpleCSsub();    % provide anything in terms of [cr, cf, cz] , get with [sr, st, sf]
        eqn = obj.getU_CS();        % provides the unit vector conversion matrix
    
            SR   you have something in terms of [sr, st, sf] and want to put it in [rx, ry, rz]
        eqn = obj.simpleSRsub();    % provide anything in terms of [sr, st, sf] , get with [rx, ry, rz]
        eqn = obj.getU_SR();        % provides the unit vector conversion matrix
        
            SC   you have something in terms of [sr, st, sf] and want to put it in [cr, cf, cz]
        eqn = obj.simpleSCsub();    % provide anything in terms of [sr, st, sf] , get with [cr, cf, cz]
        eqn = obj.getU_SC();        % provides the unit vector conversion matrix
    
        help_diff()        % ref to length, sA, and vol   for rec, cyn, sph
        help_namb()        % ref to del operator in rec, cyn, sph
        help_divg()        % ref to divergence calculation for rec, cyn, sph
        help_curl()        % ref to curl calculation for rec, cyn, sph
        help_laplacian()   % ref to scalar laplacians rec, cyn, sph    
    
        namb = getNambRec()      returns [ 1, 1, 1]        if you just want the del-operator to use
        namb = getNambCyn()      returns [ 1, 1/cr, 1]
        namb = getNambSph()      returns [ 1, 1/sr, 1/(sr*sin(st)) ]
        
        gradV = getGradRec(scal); rec scalar function --> gradient properly multiplied  VECTOR
        gradV = getGradCyn(scal); cyn scalar function --> gradient properly multiplied  VECTOR
        gradV = getGradSph(scal); sph scalar function --> gradient properly multiplied  VECTOR
    
        divg = getDivgRec(vecF)   input a rec vector field, get dot(namb, vecF) back SCALAR
        divg = getDivgCyn(vecF)   input a cyn vector field, get dot(namb, vecF) back SCALAR
        divg = getDivgSph(vecF)   input a sph vector field, get dot(namb, vecF) back SCALAR
    
        curlOut = getCurlRec(vec)  input rec vector field, get curl   VECTOR
        curlOut = getCurlCyn(vec)  input cyn vector field, get curl   VECTOR
        curlOut = getCurlSph(vec)  input sph vector field, get curl   VECTOR
    
        lapS = getLaplacianRecS(scal)     input rec scalar, get Laplacian SCALAR
        lapV = getLaplacianRecV(vec)      input rec vector, get Laplacian VECTOR 
        lapS = getLaplacianCynS(scal)     input cyn scalar, get Laplacian SCALAR
        lapV = getLaplacianCynV(vec)      input cyn vector, get Laplacian VECTOR
        lapS = getLaplacianSphS(scal)     input sph scalar, get Laplacian SCALAR
        
        UTILITY
        seg = t01seg(start, stop)    rec ONLY...str8 lines, pt[0,1]   from start to stop
        dist = distRec(pa,pb)    the distance between 2 rec points
        dist = distCyn(pa,pb)    the distance between 2 cyn points (not on arc)
        dist = distSph(pa,pb)    the distance between 2 sph points (not on arc)
        vec_out = crossRec(vec_a, vec_b) give it 2 rec vectors and get the cp
                the 6 unit transforms...no change to vin, but mixed coords:
        vout = swapRS(vin)  convert the unit vectors <arx,ary,arz> -> <acr, acf, acz>
        
        fix...starts @ line 359
        just convert to rec first if you want these...not static, needs convt
        vec_out = crossCyn(this, vec_a, vec_b) give it 2 cyn vectors and get the curl of them
        vec_out = crossSph(this, vec_a, vec_b) give it 2 sph vectors and get the curl of them
    
        CH4
        force = CLFpt(q_poi, r_poi, q_source, r_source)   % rec point charges
    
        * add: jacobians, partials, ect      use systems as bridge to eachother if needed
    %}
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
    
    properties
        pts;      % once fed, tell it input format, other 2 coordinates populate
        inType;   % stores string specifiny system input coordinates are provided in
    end
    
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
    
    methods
        
        function this = cls_PTconvert(str) % too many input arguments ?
            if nargin == 1
                printf('%s  ', str);
            end
            this.pts = nan(3, 3);
            %fprintf('\npoint converting obj made\n');
            %disp(this.pts);
        end
%------------------------------------------------------------------------------------------        
        
        function print(this)
            if isnan(this.pts(1,1))
                fprintf('it is empty\n');
            else
                fprintf('\n    *** angles are actually held in radians ***...original: %s\n',...
                    this.inType);
                fprintf('rec (x, y, z)   :  %.4f  ,  %.4f   ,  %.4f  \n', this.pts(1, :) );
                fprintf('cyn ( r, fi, z) :  %.4f  ,  %.4f°  ,  %.4f  \n',...
                    this.pts(2, 1), rad2deg(this.pts(2, 2)), this.pts(2, 3) );
                fprintf('sph (r, th, fi) :  %.4f  ,  %.4f°  ,  %.4f° \n',...
                    this.pts(3, 1), rad2deg(this.pts(3, 2)), rad2deg(this.pts(3, 3)) );
            end
        end
%------------------------------------------------------------------------------------------        
        
        function mag(this)
            fprintf('\n    *** all the magnitudes should be the same *** \n');
            fprintf('\nrec length :  %.4f \n',...
                sqrt( sum(this.pts(1, :).^2) ) ); % the clean way
            fprintf('cyn length :  %.4f \n',...
                sqrt( this.pts(2, 1)^2  + this.pts(2, 3)^2 ));
            fprintf('sph length :  %.4f \n', this.pts(3, 1) );
        end
%------------------------------------------------------------------------------------------        
        
        function clear(this)
            this.pts = nan(3, 3);
            this.inType = 'none';
        end
%------------------------------------------------------------------------------------------        
        
        function feed(this,cords, type)
            this.clear();
            if  type == 'R'
                this.inType = 'rec';
                this.pts(1, 1) = cords(1);
                this.pts(1, 2) = cords(2);
                this.pts(1, 3) = cords(3);
                this.rec_to_cyn();
                this.rec_to_sph();
            end
 
            if  type == 'C'
                this.inType = 'cyn';
                this.pts(2, 1) = cords(1);
                this.pts(2, 2) = cords(2);
                this.pts(2, 3) = cords(3);
                this.cyn_to_rec();
                this.rec_to_sph();  % ok to call after others, using rec as a bridge
            end
            
            if  type == 'S'
                this.inType = 'sph';
                this.pts(3, 1) = cords(1);
                this.pts(3, 2) = cords(2);
                this.pts(3, 3) = cords(3);
                this.sph_to_rec();
                this.rec_to_cyn(); % circular call
            end 
        end
%------------------------------------------------------------------------------------------ 
        
        function rec_to_cyn(this)
            Prx = this.pts(1, 1);
            Pry = this.pts(1, 2);
            Prz = this.pts(1, 3);
            [Pcf, Pcr, Pcz] = cart2pol(Prx, Pry, Prz);
            this.pts(2, 3) = Pcz;
            this.pts(2, 1) = Pcr;
            if Pcf >= 0
                this.pts(2, 2) = Pcf;
            else
                this.pts(2, 2) = 2*pi + Pcf;
            end
        end
%------------------------------------------------------------------------------------------        
        
        function rec_to_sph(this)
            Prx = this.pts(1, 1);
            Pry = this.pts(1, 2);
            Prz = this.pts(1, 3);
            [Psf, Pst, Psr] = cart2sph(Prx, Pry, Prz);
            if Psf >= 0
                this.pts(3, 3) = Psf;
            else
                this.pts(3, 3) = Psf + 2*pi;
            end
            if Pst >= 0
                this.pts(3, 2) = (pi/2) - Pst;
            else
                this.pts(3, 2) = (pi/2) - Pst;
            end
            this.pts(3, 1) = Psr;    
        end
%------------------------------------------------------------------------------------------        
        
        function cyn_to_rec(this)
            Pcr = this.pts(2, 1);
            Pcf = this.pts(2, 2);
            Pcz = this.pts(2, 3);
            if Pcf >= pi
                Pcf = Pcf - 2*pi;
            end
            [Prx, Pry, Prz] = pol2cart(Pcf, Pcr, Pcz);
            this.pts(1, 1) = Prx;
            this.pts(1, 2) = Pry;
            this.pts(1, 3) = Prz;
        end
%------------------------------------------------------------------------------------------        
        
        function sph_to_rec(this)
            Psr = this.pts(3, 1);
            Pst = this.pts(3, 2);
            Psf = this.pts(3, 3);
            if Psf >= pi
                Psf = Psf - 2*pi;
            end
            if Pst <= pi/2
                Pst = (pi/2) - Pst;
            else
                Pst = (pi/2) - Pst;
            end
            [Prx, Pry, Prz] = sph2cart(Psf, Pst, Psr);
            this.pts(1, 1) = Prx;
            this.pts(1, 2) = Pry;
            this.pts(1, 3) = Prz; 
        end
%------------------------------------------------------------------------------------------        
        
        function graph(this)
            figure();
            maxP = this.pts(3,1); % just take the radius
            hold on;
            axis equal;
            grid on;
            view(125, 20);
            title('CURRENT POINT', 'fontsize', 16);
            xlabel('x axis');
            ylabel('y axis');
            zlabel('z axis');
            xlim([-maxP-1, maxP+1]);
            ylim([-maxP-1, maxP+1]);
            zlim([-maxP-1, maxP+1]);
            xax = linspace(-maxP-1, maxP+1, 128);
            yax = linspace(-maxP-1, maxP+1, 128);
            zax = linspace(-maxP-1, maxP+1, 128);
            plot3(xax, 0*xax, 0*xax, 'k', 'linewidth', 1);
            plot3(0*yax, yax, 0*yax, 'k', 'linewidth', 1);
            plot3(0*zax, 0*zax, zax, 'k', 'linewidth', 1);
            plot3(maxP+1,0,0,'yo', 'markersize', 6, 'linewidth', 10); % +x direction
            text(maxP+1,0,0,'+X', 'FontSize', 12);
            plot3(0,maxP+1,0,'yo', 'markersize', 6, 'linewidth', 10); % +y direction
            text(0,maxP+1,0,'+Y', 'FontSize', 12);
            plot3(0,0,maxP+1,'yo', 'markersize', 6, 'linewidth', 10); % +z direction
            text(0,0,maxP+1,'+Z', 'FontSize', 12);
            plot3(this.pts(1, 1), this.pts(1, 2), this.pts(1, 3), 'b.', 'markersize', 20);
            plot3([0, this.pts(1, 1)], [0, this.pts(1, 2)], [0, this.pts(1, 3)],...
                'r:', 'linewidth', 2);
            strR = sprintf("rec ( %.2f, %.2f, %.2f )", this.pts(1,:));
            strC = sprintf("cyn ( %.2f, %.2f°, %.2f )",...
                this.pts(2, 1), rad2deg(this.pts(2, 2)), this.pts(2, 3) );
            strS = sprintf("sph ( %.2f, %.2f°, %.2f° )",...
                this.pts(3,1), rad2deg(this.pts(3, 2)), rad2deg(this.pts(3, 3)) );
            text(maxP,0, maxP, strR,'fontweight','bold','fontsize',13);
            text(maxP,0, maxP-1, strC,'fontweight','bold','fontsize',13);
            text(maxP,0, maxP-2, strS,'fontweight','bold','fontsize',13);
            r = this.pts(2,1);
            ht = maxP+1;
            [cordX, cordY, cordZ] = cylinder(r);
            cordZ = cordZ * ht;
            surf(cordX, cordY, cordZ,'FaceColor', 'b', 'FaceAlpha', .1, 'EdgeColor', 'none');
            cordZ = (-1) .* cordZ;
            surf(cordX, cordY, cordZ,'FaceColor', 'b', 'FaceAlpha', .1, 'EdgeColor', 'none');
            [sX, sY, sZ] = sphere;
            r = this.pts(3,1);
            surf(r.*sX, r.*sY, r.*sZ, 'FaceColor', 'g', 'FaceAlpha', .1, 'EdgeColor', 'none');   
        end 
        
%------------------------------------------------------------------------------------------
        function trans = transVecRC(this, vecIn)    % transform vec [rx, ry, rz] --> [cr, cf, cz]
            tempU = this.getU_RC();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleRCsub(tempV);         % completley transform variables
                %temp = this.simpleSCsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
        function trans = transVecRS(this, vecIn)    % transform vec [rx, ry, rz] --> [sr, st, sf]
            tempU = this.getU_RS();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleRSsub(tempV);         % completley transform variables
                %temp = this.simpleCSsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
        function trans = transVecCR(this, vecIn)    % transform vec [cr, cf, cz] --> [rx, ry, rz]
            tempU = this.getU_CR();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleCRsub(tempV);         % completley transform variables
                %temp = this.simpleSRsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
        function trans = transVecCS(this, vecIn)    % transform vec [cr, cf, cz] --> [sr, st, sf]
            tempU = this.getU_CS();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleCSsub(tempV);         % completley transform variables
                %temp = this.simpleRSsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
        function trans = transVecSR(this, vecIn)    % transform vec [sr, st, sf] --> [rx, ry, rz]
            tempU = this.getU_SR();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleSRsub(tempV);         % completley transform variables
                %temp = this.simpleCRsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
        function trans = transVecSC(this, vecIn)    % transform vec [sr, st, sf] --> [cr, cf, cz]
            tempU = this.getU_SC();                 % get the appropriate trans vector
            tempV = vecIn * transpose(tempU);       % perform vector multiplication
            temp = this.simpleSCsub(tempV);         % completley transform variables
                %temp = this.simpleRCsub(tempV);
            trans = simplify(temp);                 % std format...simplify is nice 
        end
%------------------------------------------------------------------------------------------
function vec_out = crossCyn(this, vec_a, vec_b)
    this.feed(vec_a, 'C');
    va = this.pts(1,:);
    this.feed(vec_b, 'C');
    vb = this.pts(1,:);
    v1 = cross(va,vb);
    this.feed(v1, 'R');
    vec_out = this.pts(2,:);
end
%------------------------------------------------------------------------------------------
function vec_out = crossSph(this, vec_a, vec_b)
    this.feed(vec_a, 'S');
    va = this.pts(1,:);
    this.feed(vec_b, 'S');
    vb = this.pts(1,:);
    v1 = cross(va,vb);
    this.feed(v1, 'R');
    vec_out = this.pts(3,:);
end
%------------------------------------------------------------------------------------------

    end % end methods
    
%******************************************************************************************   
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------ 
%******************************************************************************************


    methods(Static) % for utility

        
%------------------------------------------------------------------------------------------ 
%   RC   you have something in terms of [rx, ry, rz] and want to put it in [cr, cf, cz]        RC
%------------------------------------------------------------------------------------------ 
%{
        *** Var Change ***                
        rx --> cr * cos(cf)     
        ry --> cr * sin(cf)     
        rz --> cz
        *** Vec Trans Mtx ***   rec to cyn...     unit vectors in rows
                          
                Acr    cos(cf)   sin(cf)   0   Arx            ie arx = < cos(cf), sin(cf), 0 >                
                Acf =  -sin(cf)  cos(cf)   0 * Ary            multiply for comp change:
                Acz    0         0         1   Arz
                                                              Acr = Arx*cos(cf)  + Ary*sin(cf)
                                                              Acf = -Arx*sin(cf) + Ary*cos(cf)
                                                              Acz = Arz

        eqn = obj.simpleRCsub();    % provide anything in terms of [rx, ry, rz] , get with [cr, cf, cz]
        eqn = obj.getU_RC();        % provides the unit vector conversion matrix
%}
%------------------------------------------------------------------------------------------ 

        function [eqnOut] = simpleRCsub(eqnIn)
            global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            %global sr; global st; global sf; % spherical params    sr, st, sf
            %
            frx = cr * cos(cf);
            fry = cr * sin(cf);
            f1 = sqrt(rx^2 + ry^2);
            f2 = rx / sqrt(rx^2 + ry^2);
            f3 = ry / sqrt(rx^2 + ry^2);
            f4 = atan2(ry, rx); 
            eqnOut = subs(eqnIn, [ rx , ry , rz, f1, f2     , f3     , f4 ],...
                                 [ frx, fry, cz, cr, cos(cf), sin(cf), cf ]);
            %{
            frx = cr * cos(cf);
            fry = cr * sin(cf);
            frz = cz;
            eqnOut = subs(eqnIn, [rx, ry, rz], [frx, fry, frz]);
            %}
        end   
%------------------------------------------------------------------------------------------

        function unitOut = getU_RC()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
                RC_trans = [  cos(cf)  , sin(cf), 0  ;                                 
                             -sin(cf)  , cos(cf), 0  ;
                                    0  , 0      , 1 ]; 
            unitOut = RC_trans;
        end
%------------------------------------------------------------------------------------------ 
%                                       END RC
%------------------------------------------------------------------------------------------ 
%******************************************************************************************
%------------------------------------------------------------------------------------------ 
%   RS   you have something in terms of [rx, ry, rz] and want to put it in [sr, st, sf]        RS
%------------------------------------------------------------------------------------------ 
%{
        *** Var Change *** 
        rx --> sr * sin(st) * cos(sf)    
        ry --> sr * sin(st) * sin(sf)   
        rz --> sr * cos(st)
        *** Vec Trans Mtx ***   rec to sph...     unit vectors in rows

            Asr     sin(st)*cos(sf)    sin(st)*sin(sf)    cos(st)     Arx      multiply for comp change:
            Ast  =  cos(st)*cos(sf)    cos(st)*sin(sf)   -sin(st)  =  Ary
            Asf            -sin(sf)            cos(sf)          0     Arz

                Asr = Arx * sin(st) * cos(sf) + Ary * sin(st)  * sin(sf) + Arz * cos(st)
                Ast = Arx * cos(st) * cos(sf) + Ary * cos(st)  * sin(sf) - Arz * sin(st) 
                Asf = -Arx * sin(sf) + Ary * cos(sf)
                                
        eqn = obj.simpleRSsub();    % provide anything in terms of [rx, ry, rz] , get with [sr, st, sf]
        eqn = obj.getU_RS();        % provides the unit vector conversion matrix

%}
%------------------------------------------------------------------------------------------

        function [eqnOut] = simpleRSsub(eqnIn)
            global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            %
            frx = sr * sin(st) * cos(sf);
            fry = sr * sin(st) * sin(sf);
            frz = sr * cos(st);
            f1 = sqrt(rx^2 + ry^2 + rz^2);
            f2 = sqrt(rx^2 + ry^2) / sqrt(rx^2 + ry^2 + rz^2);
            f3 = rz / sqrt(rx^2 + ry^2 + rz^2);
            f4 = ry / sqrt(rx^2 + ry^2);
            f5 = rx / sqrt(rx^2 + ry^2);
            f6 = acos( rz / sqrt( rx^2 + ry^2 + rz^2 ) );
            f7 = atan2(ry, rx);
            eqnOut = subs(eqnIn, [ rx , ry , rz , f1, f2     , f3     , f4     , f5     , f6, f7 ],...
                                 [ frx, fry, frz, sr, sin(st), cos(st), sin(sf), cos(sf), st, sf ]);
            %{
            frx = sr * sin(st) * cos(sf);
            fry = sr * sin(st) * sin(sf);
            frz = sr * cos(st);
            eqnOut = subs(eqnIn, [rx, ry, rz], [frx, fry, frz]);
            %}
        end
%------------------------------------------------------------------------------------------

        function unitOut = getU_RS()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf        
            RS_trans = [  sin(st) * cos(sf) , sin(st) * sin(sf), cos(st)  ;        
                          cos(st) * cos(sf) , cos(st) * sin(sf), -sin(st) ;
                                   -sin(sf) ,           cos(sf), 0      ] ;
            unitOut = RS_trans;
        end
%------------------------------------------------------------------------------------------ 
%                                       END RS
%------------------------------------------------------------------------------------------ 
%******************************************************************************************
%------------------------------------------------------------------------------------------ 
%   CR   you have something in terms of [cr, cf, cz] and want to put it in [rx, ry, rz]        CR
%------------------------------------------------------------------------------------------ 
%{
        *** Var change ***                
        cr --> sqrt( rx^2 + ry^2 )        
        cf --> atan( ry / rx )  "atan2"
                        ----> sin(cf) = ry / sqrt ( rx^2 + ry^2 )
                        ----> cos(cf) = rx / sqrt ( rx^2 + ry^2 )
        cz --> rz
        *** Vec Trans Mtx ***   cyn to rec...     unit vectors in rows

                Arx    cos(cf)  -sin(cf)  0   Acr                            
                Ary =  sin(cf)  cos(cf)   0 * Acf             multiply for comp change:
                Arz    0        0         1   Acz

    Arx = Acr * cos(cf) - Acf * sin(cf)  =  Acr * rx / sqrt(rx^2 + ry^2) - Acf * ry / sqrt(rx^2 + ry^2)
    Ary = Acr * sin(cf) + Acf * cos(cf)  =  Acr * ry / sqrt(rx^2 + ry^2) + Acf * rx / sqrt(rx^2 + ry^2)
    Arz = Arz

        eqn = obj.simpleCRsub();    % provide anything in terms of [cr, cf, cz] , get with [rx, ry, rz]
        eqn = obj.getU_CR();        % provides the unit vector conversion matrix
%}
%------------------------------------------------------------------------------------------

        function [eqnOut] = simpleCRsub(eqnIn)
            global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            %
            fcfCos = rx / sqrt(rx^2 + ry^2);
            fcfSin = ry / sqrt(rx^2 + ry^2);
            fcr = sqrt(rx^2 + ry^2);
            fcf = atan(ry/ rx); % atan2(ry, rx)
            eqnOut = subs(eqnIn, [ cos(cf), sin(cf), cr , cz, cf  ],...
                                 [ fcfCos , fcfSin , fcr, rz, fcf ]);
            %{
            fcr = sqrt(rx^2 + ry^2);
            fcz = rz;
            fcf = atan2(ry, rx);
            fcfCos = rx / sqrt(rx^2 + ry^2);
            fcfSin = ry / sqrt(rx^2 + ry^2);
            eqnOut = subs(eqnIn, [cr, cz, cf cos(cf), sin(cf)],...
                [fcr, fcz, fcf, fcfCos, fcfSin]);
            %}

        end
%------------------------------------------------------------------------------------------

        function unitOut = getU_CR()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            CR_trans = [ cos(cf)   ,  -sin(cf), 0  ;
                         sin(cf)   ,  cos(cf) , 0  ;
                         0         ,  0       , 1 ];
            unitOut = CR_trans;
        end
%------------------------------------------------------------------------------------------ 
%                                       END CR
%------------------------------------------------------------------------------------------ 
%******************************************************************************************
%------------------------------------------------------------------------------------------ 
%   CS   you have something in terms of [cr, cf, cz] and want to put it in [sr, st, sf]        CS   
%------------------------------------------------------------------------------------------ 
%{
    *** Var change *** 
    cr --> sr * sin(st)
    cf --> sf
    cz --> sr * cos(st)
    *** Vec Trans Mtx ***   cyn to rec...     unit vectors in rows
            
                Asr    sin(st)  0    cos(st)    Acr                  
                Ast =  cos(st)  0   -sin(st) *  Acf            multiply for comp change
                Asf    0        1          0    Acz

        Asr = Acr * sin(st) + Acz * cos(st)
        Ast = Acr * cos(st) - Acz * sin(st)
        Asf = Acf
    
    eqn = obj.simpleCSsub();    % provide anything in terms of [cr, cf, cz] , get with [sr, st, sf]
    eqn = obj.getU_CS();        % provides the unit vector conversion matrix

%}
%------------------------------------------------------------------------------------------

        function [eqnOut] = simpleCSsub(eqnIn)                                
            global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            %
            f1 = sqrt( cr^2 + cz^2);
            fcr = sr*sin(st);
            fcz = sr*cos(st);
            eqnOut = subs(eqnIn, [f1, cr , cf, cz  ],...
                                 [sr, fcr, sf, fcz ]);
            %{
            fcr = sr*sin(st);
            fcf = sf;
            fcz = sr*cos(st);
            eqnOut = subs(eqnIn, [cr, cf, cz], [fcr, fcf, fcz] );
            %}
        end
%------------------------------------------------------------------------------------------

function unitOut = getU_CS()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf                
            CS_trans = [ sin(st), 0,  cos(st) ;
                         cos(st), 0, -sin(st) ;
                         0      , 1,       0 ];
            unitOut = CS_trans;
        end
%------------------------------------------------------------------------------------------
%                                       END CS
%------------------------------------------------------------------------------------------ 
%******************************************************************************************
%------------------------------------------------------------------------------------------ 
%   SR   you have something in terms of [sr, st, sf] and want to put it in [rx, ry, rz]        SR   
%------------------------------------------------------------------------------------------ 
%{
        *** Var change ***
        sr --> sqrt( rx^2 + ry^2 + rz^2 )
        st --> acos( rz / sqrt( rx^2 + ry^2 + rz^2 ) )
                    ---> cos(st) = rz / sqrt( rx^2 + ry^2 + rz^2 )
                    ---> sin(st) = sqrt( rx^2 + ry^2 ) / sqrt( rx^2 + ry^2 + rz^2 )
        sf --> atan (ry / rx)  "atan2"
                    ---> cos(sf) = rx / sqrt( rx^2 + ry^2 )
                    ---> sin(sf) = ry / sqrt( rx^2 + ry^2 )
        *** Vec Trans Mtx ***   cyn to rec...     unit vectors in rows
        
                Arx    sin(st) * cos(sf)    cos(st) * cos(sf)    -sin(sf)        Asr                  
                Ary =  sin(st) * sin(sf)    cos(st) * sin(sf)    cos(sf)    *    Ast     multiply for comps
                Arz    cos(st)              -sin(st)             0               Asf
        
    the component list is best done with a computer...it will make a mess

    eqn = obj.simpleSRsub();    % provide anything in terms of [sr, st, sf] , get with [rx, ry, rz]
    eqn = obj.getU_SR();        % provides the unit vector conversion matrix

%}
%------------------------------------------------------------------------------------------

        function [eqnOut] = simpleSRsub(eqnIn)                                
            global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            %
            fsr = sqrt(rx^2 + ry^2 + rz^2);
            fstSin = sqrt(rx^2 + ry^2) / sqrt(rx^2 + ry^2 + rz^2);
            fstCos = rz / sqrt(rx^2 + ry^2 + rz^2);
            fsfSin = ry / sqrt(rx^2 + ry^2);
            fsfCos = rx / sqrt(rx^2 + ry^2);
            fst = acos( rz / sqrt( rx^2 + ry^2 + rz^2 ) );
            fsf = atan2(ry, rx);
            eqnOut = subs(eqnIn, [ sr , sin(st), cos(st), sin(sf), cos(sf), st , sf  ],...
                                 [ fsr, fstSin , fstCos , fsfSin , fsfCos , fst, fsf ] );
            %{
            fsr = sqrt(rx^2 + ry^2 + rz^2);
            fst = acos( rz / sqrt( rx^2 + ry^2 + rz^2 ) );
            fstSin = sqrt(rx^2 + ry^2) / sqrt(rx^2 + ry^2 + rz^2);
            fstCos = rz / sqrt(rx^2 + ry^2 + rz^2);
            fsf = atan2(ry, rx);
            fsfSin = ry / sqrt(rx^2 + ry^2);
            fsfCos = rx / sqrt(rx^2 + ry^2);
            eqnOut = subs(eqnIn, [sr, st, sin(st), cos(st), sf, sin(sf), cos(sf) ],...
                [fsr, fst, fstSin, fstCos, fsf, fsfSin, fsfCos ] );
            %}
        end
%------------------------------------------------------------------------------------------

        function unitOut = getU_SR()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf                  
            SR_trans = [ sin(st) * cos(sf) , cos(st) * cos(sf) , -sin(sf);      
                         sin(st) * sin(sf) , cos(st) * sin(sf) , cos(sf) ;
                         cos(st)           , -sin(st)          , 0      ];
            unitOut = SR_trans;
        end
%------------------------------------------------------------------------------------------
%                                       END SR
%------------------------------------------------------------------------------------------ 
%******************************************************************************************
%------------------------------------------------------------------------------------------ 
%   SC   you have something in terms of [sr, st, sf] and want to put it in [cr, cf, cz]        SC   
%------------------------------------------------------------------------------------------ 
%{
        *** Var change ***
        sr --> sqrt( cr^2 + cz^2 )
        st --> atan( cr / cz )         check behavior with points below xy plane
        sf --> cf
        
        *** Vec Trans Mtx ***   cyn to rec...     unit vectors in rows
        
             Acr    sin(st)   cos(st)     0     Asr                  
             Acf =  0         0           1  *  Ast          multiply for comps...also messy
             Acz    cos(st)   -sin(st)    0     Asf

    eqn = obj.simpleSCsub();    % provide anything in terms of [sr, st, sf] , get with [cr, cf, cz]
    eqn = obj.getU_SC();        % provides the unit vector conversion matrix

%}
%------------------------------------------------------------------------------------------ 

        function [eqnOut] = simpleSCsub(eqnIn)                                 
            %global rx; global ry; global rz; % rectangular params  rx, ry, rz        
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf
            %{
            fsr = sqrt(cr^2 + cz^2);
            fst = atan(cr / cz); % this will cause a problem if point is below xy plane ??
            eqnOut = subs(eqnIn, [sr, st, sf], [fsr, fst, cf] );
            %}
            fsr = sqrt(cr^2 + cz^2);
            fst = atan(cr / cz); % this will cause a problem if point is below xy plane ??
            eqnOut = subs(eqnIn, [sr, st, sf], [fsr, fst, cf] );

        end
%------------------------------------------------------------------------------------------ 

        function unitOut = getU_SC()
            global rx; global ry; global rz; % rectangular params  rx, ry, rz
            global cr; global cf; global cz; % cylindrical params  cr, cf, cz
            global sr; global st; global sf; % spherical params    sr, st, sf 
            SC_trans = [ sin(st), cos(st) ,  0  ;
                         0      , 0       ,  1  ;
                         cos(st), -sin(st),  0 ];
            unitOut = SC_trans;
            
        end
%------------------------------------------------------------------------------------------
%                                       END SC
%------------------------------------------------------------------------------------------ 





%******************************************************************************************
%******************************************************************************************
%******************************************************************************************





%------------------------------------------------------------------------------------------
%       Help for reference: differentials, namblas, gradients, divergence               HELP
%------------------------------------------------------------------------------------------ 
function help_diff()  % where everything comes from
    fprintf("\ndon't mess with direction in length, limits of integration handle\n");
    fprintf("orientate surface with normal outside ...volume out \n");
    fprintf("if you were to walk around the region, it should always be on left hand\n");
    
    fprintf('\nREC [  arx, ary, arz ]  :\n');
    fprintf('\tlength:  <  d{rx}        ,  d{ry}        ,  d{rz}       >  VECTOR\n');
    fprintf('\tsA:      <  d{ry} d{rz}  ,  d{rx} d{rz}  ,  d{rx} d{ry} >  VECTOR\n');
    fprintf('vol:  d{rx} d{ry} d{rz}  SCALAR\n');
    
    fprintf('\nCYN [  acr, acf, acz ]  :\n');
    fprintf('\tlength:  <  d{cr}           ,  cr d{cf}     ,  d{cz}          >  VECTOR\n');
    fprintf('\tsA:      <  cr d{cf} d{cz}  ,  d{cr} d{cz}  ,  cr d{cr} d{cf} >  VECTOR\n');
    fprintf('vol:  cr d{cr} d{cf} d{cz}  SCALAR\n');
    
    fprintf('\nSPH  [ asr, ast, asf ]  :\n');
    fprintf('\tlength:  <  d{sr}                        ,   sr d{st}              , sr sin(st) d{sf}  >  VECTOR\n');
    fprintf('\tsA:      <  (sr)^2  sin(st) d{st} d{sf}  ,  sr sin(st) d{sr} d{sf} , sr d{sr} d{st}    >  VECTOR\n');
    fprintf('vol:  (sr)^2 sin(st) d{sr} d{st} d{sf}  SCALAR\n\n');
end
%------------------------------------------------------------------------------------------ 
function help_namb()   % the del operator ... use for gradients, divg, curl, Laplacian, ect...
    fprintf('\n fyi...gradient is    nambla * scalar-function   -->   VECTOR\n ');
    fprintf('differential operators, [ a1 , a2, a3 ] ...unit vectors may be dependent\n');
    fprintf('namb1rec =  <  d{rx}  ,  d{ry}         , d{rz}  >  VECTOR\n');
    fprintf('namb1cyn =  <  d{cr}  ,  (1/cr) d{cf}  , d{cz}  >  VECTOR\n');
    fprintf('namb1sph =  <  d{sr}  ,  (1/sr) d{st}  , 1 / ( sr * sin(st) ) d{sf}  >  VECTOR\n\n'); 
end
%------------------------------------------------------------------------------------------
function help_divg()
    fprintf('simply   dot( nambla, vector field )      out-fulx per unit volume   SCALAR:\n');
    fprintf('\n if the divergence = 0, the vector is "solenoidal"');
    fprintf(' divg ( curl (Av) ) = 0      divg of curl = 0\n');
    fprintf('REC:  d{rx} Arx  +  d{ry} Ary  +  d{rz} Arz    SCALAR\n');
    fprintf('CYN:  (1/cr) * d{cr}(cr * Acr)  +  (1/cr) * d{cf} Acf  +  d{cz} Acz  SCALAR\n');
    fprintf('SPH:  (1/sr^2) * d{sr} (sr^2)Asr  +  1/(sr * sin(st)) d{st} ( sin(st) * Ast)  +');
    fprintf('  1/(sr * sin(st)) * d{sf} Asf    SCALAR\n\n');
end
%------------------------------------------------------------------------------------------
function help_curl()
    fprintf('curl is a VECTOR ... = [0, 0, 0] then conservative/potential/irrational\n');
    fprintf('curl is just cross( nambla, vector)    rotation tendancy');
    fprintf('curl of scalar gradient = [0, 0, 0]     curl ( namb scalar) = [0, 0, 0]\n');
    fprintf('\nrec_comp1 :  d{ry} Arz - d{rz} Ary\n');
    fprintf('rec_comp2 :  d{rz} Arx  -  d{rx} Arz\n');
    fprintf('rec_comp3 :  d{rx} Ary  -  d{ry} Arx\n');
    fprintf('\ncyn_comp1  :  (1/cr) * d{cf} Acz  -  d{cz} Acf\n');
    fprintf('cyn_comp2  :  d{cz} Acr  -  d{cr} Acz\n');
    fprintf('cyn_comp3  :  (1/cr) * [ d{cr} cr * Acf  -  d{cf} Acr ]\n');
    fprintf('\nsph_comp1  :  (1/ (sr * sin(st) ) * [ d{st} sin(st) * Asf  -  d{sf} Ast ]\n');
    fprintf('sph_comp2  :  (1/sr) * [ ( 1/sin(st) ) * d{sf} Asr  -  d{sr} sr * Asf ] \n');
    fprintf('sph_comp3  :  (1/sr) * [ d{sr} sr * Ast  -  d{st} Asr ]\n');
end
%------------------------------------------------------------------------------------------
function help_laplacian()
    fprintf('\nLaplacian is a double nambla diff operator  if scalar-->scalar, vector-->vector\n');
    fprintf('a Laplacian vector is a scalar laplacian on each component');
    fprintf('if the scalar laplacian produces a value of 0, it is harmonic\n');
    fprintf('\nrec scalar:  d2{rx} V  +  d2{ry} V  +  d2{rz} V        SCALAR\n');
    fprintf('rec vector:  <  Namb2 Arx  ,  Namb2 Ary  ,  Namb2 Arz  >     VECTOR\n');
    fprintf('cyn scalar:  (1/cr) * d{cr} [ cr * d{cr} V ]  +  1/(cr^2) * d2{cf} V   +  d2{cz} V     SCALAR\n');
    fprintf('sph scalar:  (1/sr^2) * d{sr} [ (sr^2) d{sr} V ]  +\n');
                          fprintf('\t     ( 1/(sr^2 * sin(st) ) * d{st} [ sin(st) d{st} V ] +\n');
                          fprintf('\t     ( 1/(sr^2 * sin(st)^2 ) * d2{sf} V\n');
    fprintf('the vector Laplacian definitions for cyn и sph are very naughty\n');
end    
%------------------------------------------------------------------------------------------ 
%                                       END Help
%------------------------------------------------------------------------------------------





%******************************************************************************************
%******************************************************************************************
%******************************************************************************************





%------------------------------------------------------------------------------------------
% Namblas are the differential 'del' operator for gradient, divergence, curl, and Laplacian  Nams
%   the 5 nambla operations (2 Laplacian techniques) are all point funtions
%                       differentiation is always at a point
% there are also material derivatives, tensors, and other crazy shit you can do
%
% when in doubt: 
%     circulation     don't take surface integral over closed surface ==> int { vol ( divg(Av) ) }
%      out flux       don't take line integral on loop ==> int { surf ( curl(Av) ) }
%------------------------------------------------------------------------------------------  
function namb = getNambRec()
    namb = [1, 1, 1];
end
%------------------------------------------------------------------------------------------
function grad = getGradRec(scalarFun)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
        %namb = [1, 1, 1];
        %temp = gradient(scalarFun, [rx, ry, rz]);
        %grad = namb .* transpose(temp);
    g1 = diff(scalarFun, rx, 1);
    g2 = diff(scalarFun, ry, 1);
    g3 = diff(scalarFun, rz, 1);
    grad = [g1, g2, g3];
end
%------------------------------------------------------------------------------------------
function divg = getDivgRec(vecF)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    divg1 = diff(vecF(1), rx, 1);
    divg2 = diff(vecF(2), ry, 1);
    divg3 = diff(vecF(3), rz, 1);
    divg = divg1 + divg2 + divg3;
end
%------------------------------------------------------------------------------------------
function curlOut = getCurlRec(vec)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    curl1 = diff( vec(3), ry, 1) - diff( vec(2), rz, 1);
    curl2 = diff( vec(1), rz, 1) - diff( vec(3), rx, 1);
    curl3 = diff( vec(2), rx, 1) - diff( vec(1), ry, 1);
    curlOut = [ curl1, curl2, curl3 ];
end
%------------------------------------------------------------------------------------------
function lapS = getLaplacianRecS(scal)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    comp1 = diff(scal, rx, 2);
    comp2 = diff(scal, ry, 2);
    comp3 = diff(scal, rz, 2);
    lapS = comp1 + comp2 + comp3;
end
%------------------------------------------------------------------------------------------
function lapV = getLaplacianRecV(vec)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    comp1 = diff( vec(1), rx, 2);
    comp2 = diff( vec(2), ry, 2);
    comp3 = diff( vec(3), rz, 2);
    lapV = [comp1, comp2, comp3];
end
%------------------------------------------------------------------------------------------

%******************************************************************************************

%------------------------------------------------------------------------------------------ 
function namb = getNambCyn()
    global cr;
    namb = [1, 1/cr, 1];
end
%------------------------------------------------------------------------------------------
function grad = getGradCyn(scalarFun)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
        %namb = [1, 1/cr, 1];
        %temp = gradient(scalarFun, [cr, cf, cz]);
        %grad = namb .* transpose(temp);
    g1 = diff(scalarFun, cr, 1);
    g2 = (1/cr) * diff(scalarFun, cf, 1);
    g3 = diff(scalarFun, cz, 1);
    grad = [g1, g2, g3];
end
%------------------------------------------------------------------------------------------
function divg = getDivgCyn(vecF)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    divg1 = (1/cr) * diff(cr * vecF(1), cr, 1);
    divg2 = (1/cr) * diff(vecF(2), cf, 1);
    divg3 = diff(vecF(3), cz, 1);
    divg = divg1 + divg2 + divg3;
end
%------------------------------------------------------------------------------------------
function curlOut = getCurlCyn(vec)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    curl1 = (1/cr) * diff( vec(3), cf, 1) - diff( vec(2), cz, 1);
    curl2 = diff( vec(1), cz, 1) - diff( vec(3), cr, 1);
    curl3 = (1/cr) * ( diff( cr * vec(2), cr, 1) - diff( vec(1), cf, 1) );
    curlOut = [ curl1, curl2, curl3 ];
end
%------------------------------------------------------------------------------------------
function lapS = getLaplacianCynS(scal)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    comp1 = (1/cr) * diff( cr * diff( scal, cr, 1), cr, 1);
    comp2 = (1/(cr^2)) * diff( scal, cf, 2);
    comp3 = diff(scal, cz, 2);
    lapS = comp1 + comp2 + comp3;
end
%------------------------------------------------------------------------------------------
function lapV = getLaplacianCynV(vec)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    
    temp1a = (1/cr) * diff( cr * diff( vec(1), cr, 1), cr, 1);
    temp1b = (1/(cr^2)) * diff( vec(1), cf, 2);
    temp1c = diff(vec(1), cz, 2);
    temp1d = temp1a + temp1b + temp1c; % laplacian of Acr
    temp1e = vec(1) / (cr^2);          % Acr / cr^2
    temp1f = (2/cr^2) * diff(vec(2), cf, 1); % (2/cr^2) * d{cf} Acf
    comp1 = temp1d - temp1e - temp1f;
    
    temp2a = (1/cr) * diff( cr * diff( vec(2), cr, 1), cr, 1);
    temp2b = (1/(cr^2)) * diff( vec(2), cf, 2);
    temp2c = diff(vec(2), cz, 2);
    temp2d = temp2a + temp2b + temp2c; % laplacian of Acf
    temp2e = vec(2) / (cr^2);          % Acf / cr^2
    temp2f = (2/cr^2) * diff(vec(1), cf, 1); % (2/cr^2) * d{cf} Acr
    comp2 = temp2d - temp2e + temp2f;
    
    temp3a = (1/cr) * diff( cr * diff( vec(3), cr, 1), cr, 1);
    temp3b = (1/(cr^2)) * diff( vec(3), cf, 2);
    temp3c = diff(vec(3), cz, 2);
    comp3 = temp3a + temp3b + temp3c; % laplacian of Acz
    
    lapV = [comp1, comp2, comp3];
end
%------------------------------------------------------------------------------------------

%******************************************************************************************

%------------------------------------------------------------------------------------------
function namb = getNambSph()
    global sr; global st;
    namb = [1, 1/sr, 1/(sr*sin(st))];
end
%------------------------------------------------------------------------------------------ 
function grad = getGradSph(scalarFun)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
        %namb = [1, 1/sr, 1/(sr*sin(st))];
        %temp = gradient(scalarFun, [sr, st, sf]);
        %grad = namb .* transpose(temp);
    g1 = diff(scalarFun, sr, 1);
    g2 = (1/sr) * diff(scalarFun, st, 1);
    g3 = (1/(sr*sin(st))) * diff(scalarFun, sf, 1);
    grad = [g1, g2, g3];
end
%------------------------------------------------------------------------------------------
function divg = getDivgSph(vecF)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    divg1 = (1/sr^2) * diff( (sr^2) * vecF(1), sr, 1);
    divg2 = (1/(sr*sin(st))) * diff( (sin(st)) * vecF(2), st, 1);
    divg3 = (1/(sr*sin(st))) * diff(vecF(3), sf, 1);
    divg = divg1 + divg2 + divg3;
end
%------------------------------------------------------------------------------------------
function curlOut = getCurlSph(vec)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    curl1 = ( 1 / (sr * sin(st)) ) * ( diff( sin(st) * vec(3), st, 1) - diff( vec(2), sf, 1) );
    curl2 = (1/sr) * ( (1/sin(st)) * diff( vec(1), sf, 1) - diff( sr * vec(3), sr, 1) );
    curl3 = (1/sr) * ( diff( sr * vec(2), sr, 1) - diff( vec(1), st, 1) );
    curlOut = [ curl1, curl2, curl3 ];
end
%------------------------------------------------------------------------------------------
function lapS = getLaplacianSphS(scal)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    comp1 = (1 / sr^2) * diff( (sr^2) * diff(scal, sr, 1), sr, 1);
    comp2 = (1/( (sr^2)*sin(st) )) * diff( (sin(st)) * diff(scal, st, 1), st, 1);
    comp3 = (1/((sr^2)*sin(st)^2)) * diff(scal, sf, 2);
    lapS = comp1 + comp2 + comp3;
end
%------------------------------------------------------------------------------------------
function lapV = getLaplacianSphV(vec)
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    
    temp1a = (1 / sr^2) * diff( (sr^2) * diff(vec(1), sr, 1), sr, 1);
    temp1b = (1/( (sr^2)*sin(st) )) * diff( (sin(st)) * diff(vec(1), st, 1), st, 1);
    temp1c = (1/((sr^2)*sin(st)^2)) * diff(vec(1), sf, 2);
    temp1d = temp1a + temp1b + temp1c; % Laplacian of Asr , first term
    temp1e = ( 2 * vec(1) / sr^2); % second term
    temp1f = (2/(sin(st)*sr^2)) * diff(sin(st)*vec(2), st, 1); % third term
    temp1g = (2/(sin(st)*sr^2)) * diff(vec(3), sf, 1); % fourth term
    comp1 = temp1d - temp1e - temp1f - temp1g;
    
    temp2a = (1 / sr^2) * diff( (sr^2) * diff(vec(2), sr, 1), sr, 1);
    temp2b = (1/( (sr^2)*sin(st) )) * diff( (sin(st)) * diff(vec(2), st, 1), st, 1);
    temp2c = (1/((sr^2)*sin(st)^2)) * diff(vec(2), sf, 2);
    temp2d = temp2a + temp2b + temp2c; % Laplacian of Ast  , first term
    temp2e = vec(2) / ( (sr^2) * (sin(st))^2 ); % second term
    temp2f = (2/sr^2) * diff(vec(1), st, 1);  % third term
    temp2g = ((2 * cos(st)) / ((sr^2) * (sin(st)^2))) * diff(vec(3), sf, 1); % fourth term
    comp2 = temp2d - temp2e + temp2f + temp2g;
    
    temp3a = (1 / sr^2) * diff( (sr^2) * diff(vec(3), sr, 1), sr, 1);
    temp3b = (1/( (sr^2)*sin(st) )) * diff( (sin(st)) * diff(vec(3), st, 1), st, 1);
    temp3c = (1/((sr^2)*sin(st)^2)) * diff(vec(3), sf, 2);
    temp3d = temp3a + temp3b + temp3c; % Laplacian of Asf  , first term
    temp3e = vec(3) / ( (sr^2) * (sin(st))^2 ); % second term
    temp3f = (2/(sin(st)*sr^2)) * diff(vec(1), sf, 1); % third term
    temp3g = ((2 * cos(st)) / ((sr^2) * (sin(st)^2))) * diff(vec(2), sf, 1); % fourth term
    comp3 = temp3d - temp3e + temp3f + temp3g;
    
    lapV = [comp1, comp2, comp3];
end
%------------------------------------------------------------------------------------------ 
%                                       END Nams
%------------------------------------------------------------------------------------------





%******************************************************************************************
%******************************************************************************************
%******************************************************************************************





%{
    general utility , operations that show up a lot

%}
%------------------------------------------------------------------------------------------
function seg = t01seg(start, stop)
    %{
        pt[0,1] given a start and stop as const rec, parameterized seg returned
            ideally:
                        Ev = [ f(rx, ry, rz) , f(rx, ry, rz) , f(rx, ry, rz) ]
                        path = ee.t01seg(start_pt, stop_pt);
                        path_d = diff(path, pt, 1);
                        temp = subs(Ev, [rx, ry, rz] , path);
                        intg = dot( temp, path_d);
                        answ = int(intg, pt, 0, 1);
    %}
    global pt;
    seg = ((1-pt) .* start) + (pt .* stop);
end
%------------------------------------------------------------------------------------------
function dist = distRec(pa,pb)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    dist = sqrt( (pa(1)-pb(1))^2 + (pa(2)-pb(2))^2 + (pa(3)-pb(3))^2 );
end
%------------------------------------------------------------------------------------------
function dist = distCyn(pa,pb)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    dist = sqrt( pa(1)^2 + pb(1)^2 - 2 * pa(1) * pb(1) * cos( pa(2)-pb(2) ) + ( pa(3) - pb(3) )^2 );
end
%------------------------------------------------------------------------------------------
function dist = distSph(pa,pb)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    d1 = pa(1)^2 + pb(1)^2;
    d2 = -2 * pa(1) * pb(1);
    d3 = sin(pa(2)) * sin(pb(2)) * cos(pa(3)-pb(3)) + cos(pa(2)) * cos(pb(2));
    dist = sqrt(  d1 + d2 * d3 );
end
%------------------------------------------------------------------------------------------
function vec_out = crossRec(vec_a, vec_b)
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
    comp_1 = ( vec_a(2) * vec_b(3) ) - ( vec_a(3) * vec_b(2) ); % arx component
    comp_2 = ( vec_a(3) * vec_b(1) ) - ( vec_a(1) * vec_b(3) ); % ary component
    comp_3 = ( vec_a(1) * vec_b(2) ) - ( vec_a(2) * vec_b(1) ); % arz component
    vec_out = [ comp_1 , comp_2 , comp_3 ];
end
%------------------------------------------------------------------------------------------
function vout = swapRS(vin)  %  [ Arx arx, Ary ary, Arz arz] -> [ Asr asr, Ast ast, Asf, asf]
    global rx; global ry; global rz;        
    global cr; global cf; global cz;  
    global sr; global st; global sf;
    global arx; global ary; global arz; 
    global acr; global acf; global acz; 
    global asr; global ast; global asf;
    global Arx; global Ary; global Arz;
    global Acr; global Acf; global Acz;
    global Asr; global Ast; global Asf;
    
    
    
end
%------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------ 
%                                       END general utility
%------------------------------------------------------------------------------------------



%******************************************************************************************
%******************************************************************************************
%******************************************************************************************


%{
    Chapter 4: Electrostatic Fields         won't have everything...know fundamentals
        *** in general directional vector is  [ poi - source ]
        coulomb's force    F = Q1 Q2 / 4 pi ep0 R^2     aR       point charges
        E as force / unit vol -->    E = Q / 4 pi ep0 R^2  aR
        line:    Q = intL( rho_L dl )    
            oo line:  E = rho_L / 2 pi ep0 rho    acr
        surface: Q = intS( rho_S dS )    
            oo sheet:          E = rho_S / 2 ep0    an
            oo parallel plate: E = rho_S / ep0    an
        volume:  Q = intV( rho_V dv )

        Gauss flux = Q_enc = intS( dot(D, dS) ) = intV( rho_v dv )      Dn = Q_enc / S
        
        maxwell_static_1   rho_v = divg(D)        D = ep0 E
        maxwell_static_2   curl(E) = [0,0,0]      E = - grad(V)
%}
%------------------------------------------------------------------------------------------
function force = CLFpt(q_poi, r_poi, q_source, r_source)
    const = cls_CONST();  % an obj instantiated in another...it works
    % force that "source" puts on point of interest "poi"
    R = r_poi - r_source; % direction source to poi
    mult = (q_poi * q_source) / ( 4 * sym(pi) * const.ep0 * norm(R)^3);
    force = mult .* R;
end



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

        
    end % end static methods
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------    
end % end class

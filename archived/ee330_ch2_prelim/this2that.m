function [ outV, inp2out, partials, jacb] = this2that(vecIn, type)
% uses jacobian to get 3x3 identity RC, RS, CR, CS, SC, SR         vector trans ...super predator
    
    flag = 1;               % if bad input
    
    syms rx;                % rectangular params  rx, ry, rz
    assume(rx, 'real');
    syms ry;
    assume(ry,'real');
    syms rz;
    assume(rz, 'real');
    syms cr;                % cylindrical params  cr, cf, cz
    assume(cr, 'real');
    syms cf;
    assume(cf,'real');
    syms cz;
    assume(cz, 'real');
    syms sr;                % spherical params  sr, st, sf
    assume(sr, 'real');
    syms st;
    assume(st,'real');
    syms sf;
    assume(sf, 'real');
    
    
    if type == double('RC')     % rec to cyn needed
        
        flag = 0;
    end
    
    
    if type == double('RS')     % rec to sph needed  
        flag = 0;
    end
    
    
    if type == double('CR')     % cyn to rec 
%{
    *** Var change ***
rx = cr * cos(cf)     
ry = cp * sin(cf)     
rz = cz             

cr = sqrt( rx^2 + ry^2)        
cf = atan( ry / rx )  
                ----> sin(cf) = ry / sqrt ( rx^2 + ry^2 )
                ----> cos(cf) = rx / sqrt ( rx^2 + ry^2 )
cz = rz
        
    *** Comp change ***
Acr = Arx * cos(cf) + Ary * sin(cf)    
Acf = -Arx * sin(cf) * Ary * cos(cf)  
Acz = Arz    
        
Arx = {Acr * rx / sqrt( rx^2 + ry^2 )} - {Acf * ry / sqrt( rx^2 + ry^2 )}                         
Ary = { Acr * ry / sqrt( rx^2 + ry^2 ) } + { Acf * rx / sqrt( rx^2 + ry^2 )} 
Arz = Acz
%}        
        frx = cr * cos(cf);
        partials(1, 1) = diff(frx, cr, 1); 
        partials(1, 2) = diff(frx, cf, 1);
        partials(1, 3) = diff(frx, cz, 1);
        fry = cr * sin(cf);
        partials(2, 1) = diff(fry, cr, 1);
        partials(2, 2) = diff(fry, cf, 1);
        partials(2, 3) = diff(fry, cz, 1);
        frz = cz;
        partials(3, 1) = diff(frz, cr, 1);
        partials(3, 2) = diff(frz, cf, 1);
        partials(3, 3) = diff(frz, cz, 1);
        
        jacb = simplify(abs(det(partials))); %  don't forget d_cr d_cf d_cz 
        inp2out = [ cos(cf), -sin(cf), 0  ;         % cyn vector -> rec vector
                    sin(cf),  cos(cf), 0  ;
                    0      ,  0      , 1 ];
        out = inp2out * transpose(vecIn);
        display(out); % raw, not simplified
        
        outV = subs( out, [ cz, cr, sin(cf), cos(cf) ],...
        [rz, sqrt( rx^2 + ry^2 ), ry / sqrt( rx^2 + ry^2 ) , rx / sqrt( rx^2 + ry^2 ) ] );
        
        flag = 0;
    end
    
    
    if type == double('CS')     % cyn to sph needed       fyi  st = atan(cr / cz)
        flag = 0;
    end
    
    
    if type == double('SR')     % sph to rec 
%{
   *** Var change ***
rx = sr * sin(st) * cos(sf)    
ry = sr * sin(st) * sin(sf)   
rz = sr * cos(st)                             
                                   
sr = sqrt( rx^2 + ry^2 + rz^2)       
st = acos( rz / sqrt( rx^2 + ry^2 + rz^2) ) 
                         ----> cos(st) = rz / sqrt( rx^2 + ry^2 + rz^2)
                         ----> sin(st) = sqrt( rx^2 + ry^2 ) / sqrt( rx^2 + ry^2 + rz^2)
sf = atan( ry / rx )
            ----> cos(sf) = rx / sqrt ( rx^2 + ry^2 )
            ----> sin(sf) = ry / sqrt ( rx^2 + ry^2 )
        
        
  *** Comp change ***
Asr = Arx * sin(st) * cos(sf) + Ary * sin(st)  * sin(sf) + Arz * cos(st)
Ast = Arx * cos(st) * cos(sf) + Ary * cos(st)  * sin(sf) - Arz * sin(st) 
                           
Arx = {Acr * rx / sqrt( rx^2 + ry^2 )} - {Acf * ry / sqrt( rx^2 + ry^2 )}
Arz = Acz
Ary = { Acr * ry / sqrt( rx^2 + ry^2 ) } + { Acf * rx / sqrt( rx^2 + ry^2 )}
%}
        frx = sr * sin(st) * cos(sf);
        partials(1, 1) = diff(frx, sr, 1); 
        partials(1, 2) = diff(frx, st, 1); 
        partials(1, 3) = diff(frx, sf, 1);
        fry = sr * sin(st) * sin(sf);
        partials(2, 1) = diff(fry, sr, 1);
        partials(2, 2) = diff(fry, st, 1);
        partials(2, 3) = diff(fry, sf, 1);
        frz = sr * cos(st);
        partials(3, 1) = diff(frz, sr, 1);
        partials(3, 2) = diff(frz, st, 1);
        partials(3, 3) = diff(frz, sf, 1);
        
        jacb = simplify(abs(det(partials)));   %  don't forget d_sr d_st d_sf   
        inp2out = [ sin(st) * cos(sf) , cos(st) * cos(sf) , -sin(sf);   % sph vec to rec vec
                    sin(st) * sin(sf) , cos(st) * sin(sf) , cos(sf) ;
                    cos(st)           , -sin(st)          , 0      ];
        
        out = inp2out * transpose(vecIn);
        display(out); % raw, not simplified
        
        outV = subs( out, [ sr, cos(st), sin(st), cos(sf), sin(sf) ],...
        [ sqrt( rx^2 + ry^2 + rz^2),...
          rz / sqrt( rx^2 + ry^2 + rz^2),...
          sqrt( rx^2 + ry^2 ) / sqrt( rx^2 + ry^2 + rz^2),...
          rx / sqrt( rx^2 + ry^2 ), ry  / sqrt( rx^2 + ry^2 ) ] );
          
        flag = 0;
    end
    
    
    if type == double('SC')     % sph to cyn needed
        flag = 0;
    end
    
    if flag == 1
        fprintf('there was a problem with type...use RC, RS, CR, CS, SR, SC only\n');
        inp2out = zeros(3, 3);
        outV = zeros(3, 1)
    end
end












%{


First to work:

function [c2r] = cyn2recV(cynd)
    c2r = zeros(3, 1);
    syms rx;                % rectangular params  rx, ry, rz
    assume(rx, 'real');
    syms ry;
    assume(ry,'real');
    syms rz;
    assume(rz, 'real');              
    syms cr;                % cylindrical params  cr, cf, cz
    assume(cr, 'real');
    syms cf;
    assume(cf,'real');
    syms cz;
    assume(cz, 'real');
    
    
    cynd = transpose(cynd);
    cynTOrec = [ cos(cf), -sin(cf), 0  ;         % cyn vector -> rec vector
                 sin(cf),  cos(cf), 0  ;
                 0      ,  0      , 1 ];
    c2r = cynTOrec * cynd;
    %display(c2r);
    
    %temp = subs(cynTOrec(1, 1), [ cr, cf, cz, rx, ry, rz, sin(cf), cos(cf) ] ,...
     %   [ cr, cf, rz, rx, ry       rz, ry / sqrt( rx^2 + ry^2 ) , rx / sqrt( rx^2 + ry^2 ) ] );
    
    temp = subs(c2r(1, 1), [ cz, cr, sin(cf), cos(cf) ],...
        [rz, sqrt( rx^2 + ry^2 ), ry / sqrt( rx^2 + ry^2 ) , rx / sqrt( rx^2 + ry^2 ) ] );
    
    Ax_eq = temp;
    
    
    display(Ax_eq);
    c2r = c2r + 0;
end



    Ax_eq = subs(cynTOrec(1, 1), [ cr, cf, cz ],...
        [ sqrt( rx^2 + ry^2 ), atan( ry /rx), rz ] );
 %Ax_eq = subs(cynTOrec(1, 1), [ rx, ry, rz, cr, cf, cz ],...
     %   [ rx, ry, rz, sqrt( rx^2 + ry^2 ), atan( ry /rx), rz ] );
%}


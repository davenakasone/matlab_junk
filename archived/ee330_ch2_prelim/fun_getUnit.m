function [unitOut] = fun_getUnit(type)
%{

        BEST VERSION...but try to use it as a static method in cls_PTconvert

 you have a vector field in R, C, or S
 you need a UNIT vector interms of R, C, or S in terms of the original vector field
 that way you can multiple the UNIT vector field output and proceed with calculation
 set 'type' to 'RC', 'RS', 'CR', 'CS', 'SR', or 'SC' to get the unit vector you desire
 unitOut will be a (3 x 3) with each component in terms of original variable

    ie, you have a vector field in terms of CYN, say Avec = [ cr^2 , sin(cf) * cos(cf), 3 * cr ]
            you want to this field to correspond to REC, but in CYN params

        conv = fun_getUnits('CR')   % conv assigned 'unitOut'  a (3x3) in terms of (cr, cf, cz)
            to represent this vector field in REC
                field = conv .* Avec;

                  unitOut will look like:   arx:  f( cr, cf, cz) , f( cr, cf, cz) , f( cr, cf, cz)
                                            ary:  f( cr, cf, cz) , f( cr, cf, cz) , f( cr, cf, cz)
                                            arz:  f( cr, cf, cz) , f( cr, cf, cz) , f( cr, cf, cz)
 if you need arx, just take first row...
%}   
    %global rx; global ry; global rz; % rectangular params  rx, ry, rz
    %global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    %global sr; global st; global sf; % spherical params    sr, st, sf
    %sym rx;  sym ry;  sym rz; 
    %sym cr;  sym cf;  sym cz; 
    %sym sr;  sym st;  sym sf; 
    
    flag = 1; % for bad input
    global rx; global ry; global rz; % rectangular params  rx, ry, rz
    global cr; global cf; global cz; % cylindrical params  cr, cf, cz
    global sr; global st; global sf; % spherical params    sr, st, sf
                        
    if type == double('RC')
        %{ 
                                                            rec to cyn     RC
            acr    cos(cf)   sin(cf)   0   arx                     
            acf =  -sin(cf)  cos(cf)   0 * ary
            acz    0         0         1   arz
        %}
            RC_trans = [ cos(cf)  , sin(cf), 0 ;                                 
                        -sin(cf) , cos(cf), 0  ;
                        0       , 0      , 1 ]; 
        unitOut = RC_trans;
        flag = 0;
    end
    
    if type == double('RS')
        %{ 
                                                            rec to sph     RS
            asr    cos(cf)  -sin(cf)  0   arx                     
            ast =  sin(cf)  cos(cf)   0 * ary
            asf    0        0         1   arz
        %}         
        RS_trans = [ sin(st) * cos(sf) , sin(st) * sin(sf), cos(st)  ;        
                    cos(st) * cos(sf)  , cos(st) * sin(sf), -sin(st) ;
                    -sin(sf)           , cos(sf)          ,  0     ] ;
        unitOut = RS_trans;
        flag = 0;
    end
    
    if type == double('CR')
        %{ 
                                                            cyn to rec     CR
            arx    cos(cf)  -sin(cf)  0   acr                     
            ary =  sin(cf)  cos(cf)   0 * acf
            arz    0        0         1   acz
        %}
        CR_trans = [ cos(cf)   ,  -sin(cf), 0  ;
                     sin(cf)   ,  cos(cf) , 0  ;
                     0         ,  0       , 1 ];
        unitOut = CR_trans;
        flag = 0;
    end
    
    if type == double('CS')
        %{
                                                            cyn to sph     CS
            asr    sin(st)  0    cos(st)    acr                  
            ast =  cos(st)  0   -sin(st) *  acf
            asf    0        1          0    acz
        %}
        CS_trans = [ sin(st), 0,  cos(st) ;
                     cos(st), 0, -sin(st) ;
                     0      , 1,       0 ];
        unitOut = CS_trans;
        flag = 0;
    end
    
    if type == double('SR')
        %{
                                                            sph to rec     SR
            arx    sin(st) * cos(st)    cos(st) * cos(sf)    -sin(sf)        asr                  
            ary =  sin(st) * sin(sf)    cos(st) * sin(sf)    cos(sf)    *    ast
            arz    cos(st)              -sin(st)             0               asf
        %} 
        SR_trans = [ sin(st) * cos(sf) , cos(st) * cos(sf) , -sin(sf);      
                     sin(st) * sin(sf) , cos(st) * sin(sf) , cos(sf) ;
                     cos(st)           , -sin(st)          , 0      ];
        unitOut = SR_trans;
        flag = 0;
    end
    
    if type == double('SC')
        %{
                                                            sph to cyn     SC
            acr    sin(st)  0    cos(st)    asr                  
            acf =  cos(st)  0   -sin(st) *  ast
            acz    0        1          0    asf
        %}
        SC_trans = [ sin(st), 0, cos(st)  ;
                     cos(st), 0, -sin(st) ;
                     0      , 1,       0 ];
        unitOut = SC_trans;
        flag = 0;
    end
    
    if flag == 1
        fprintf('\n     the input type was not properly specified\n');
        unitOut = nan(3, 3);
    end
    %display(unitOut); % turn off or on here
end




function [vecs] = fun_vecCon(pts, symVec, type)
% get that 3x3 with parameters in 'points' matrix   book format
% symVec as definition of one of the vectors
% type R, C, S   match to vector

    vecs = zeros(3, 3);     % keep first 3x3 for symbols, last 3x3 for evaluating points
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
    
    recTOcyn = [ cos(cf) , sin(cf), 0;                                   % RCtrans
                 -sin(cf), cos(cf), 0;
                 0       , 0      , 1];
    
    recTOsph = [ sin(st) * cos(sf) , sin(st) * sin(sf), cos(st) ;        % RStrans
                 cos(st) * cos(sf) , cos(st) * sin(sf), -sin(st);
                 -sin(sf)          , cos(sf)          ,  0     ];
             
    cynTOrec = [ cos(cf), -sin(cf), 0  ;                                 % CRtrans
                 sin(cf),  cos(cf), 0  ;
                 0      ,  0      , 1 ];
    
             
             
             % need CS
             
             
    sphTOrec = [ sin(st) * cos(sf) , cos(st) * cos(sf) , -sin(sf);       % SRtrans
                 sin(st) * sin(sf) , cos(st) * sin(sf) , cos(sf) ;
                 cos(st)           , -sin(st)          , 0      ];
          
    sphTOcyn = [ sin(st) , cos(st) , 0;                                  %SCtrans
                 0       , 0       , 1;
                 cos(st) , -sin(st), 0];
    
    if type == 'R'
        Arx = symVec(1, 1);
        Ary = symVec(1, 2);
        Arz = symVec(1, 3);
        ArS = [ Arx; Ary; Arz ];
        vecs(1, 1) = subs(Arx, [rx, ry, rz] , [pts(1, 1), pts(1, 2), pts(1, 3)]);
        vecs(1, 2) = subs(Ary, [rx, ry, rz] , [pts(1, 1), pts(1, 2), pts(1, 3)]);
        vecs(1, 3) = subs(Arz, [rx, ry, rz] , [pts(1, 1), pts(1, 2), pts(1, 3)]);
        
        RCtrans = recTOcyn * ArS; % functions of rx, ry, rz, cr, cf, cz
        fprintf('\nA_cr : %s \n', RCtrans(1, 1));
        vecs(2, 1) = subs(RCtrans(1, 1), [ rx, ry, rz, cr, cf, cz ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(2, 1), pts(2, 2), pts(2, 3) ] );
        fprintf('A_cf : %s \n', RCtrans(2, 1));
        vecs(2, 2) = subs(RCtrans(2, 1), [ rx, ry, rz, cr, cf, cz ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(2, 1), pts(2, 2), pts(2, 3) ] );
        fprintf('A_cz : %s \n', RCtrans(3, 1));
        vecs(2, 3) = subs(RCtrans(3, 1), [ rx, ry, rz, cr, cf, cz ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(2, 1), pts(2, 2), pts(2, 3) ] );
        
        RStrans = recTOsph * ArS; % functions of rx, ry, rz, sr, st, sf
        fprintf('\nA_sr : %s \n', RStrans(1, 1));
        vecs(3, 1) = subs(RStrans(1, 1), [ rx, ry, rz, sr, st, sf ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(3, 1), pts(3, 2), pts(3, 3) ] );
        fprintf('A_st : %s \n', RStrans(2, 1));
        vecs(3, 2) = subs(RStrans(2, 1), [ rx, ry, rz, sr, st, sf ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(3, 1), pts(3, 2), pts(3, 3) ] );
        fprintf('A_sf : %s \n', RStrans(3, 1));
        vecs(3, 3) = subs(RStrans(3, 1), [ rx, ry, rz, sr, st, sf ],...
            [ pts(1, 1), pts(1, 2), pts(1, 3), pts(3, 1), pts(3, 2), pts(3, 3) ] );
    end
    
    if type == 'C'
        Acr = symVec(1, 1);
        Acf = symVec(1, 2);
        Acz = symVec(1, 3);
        AcS = [Acr; Acf; Acz];
        
        
        
        
        CRtrans = cynTOrec * AcS; % functions of cr, cf, cz, rx, ry, rz
        
        fprintf('\nA_cr : %s \n', CRtrans(1, 1));
        vecs(1, 1) = subs(CRtrans(1, 1), [ cr, cf, cz, rx, ry, rz ],...
            [ pts(2, 1), pts(2, 2), pts(2, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_cf : %s \n', CRtrans(2, 1));
        vecs(1, 2) = subs(CRtrans(2, 1), [ cr, cf, cz, rx, ry, rz ],...
            [ pts(2, 1), pts(2, 2), pts(2, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_cz : %s \n', CRtrans(3, 1));
        vecs(1, 3) = subs(CRtrans(3, 1), [ cr, cf, cz, rx, ry, rz ],...
            [ pts(2, 1), pts(2, 2), pts(2, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        
        
    end
    
    
    
    if type == 'S'
        Asr = symVec(1, 1);
        Ast = symVec(1, 2);
        Asf = symVec(1, 3);
        AsS = [Asr; Ast; Asf]; 
        vecs(3, 1) = subs(Asr, [sr, st, sf] , [pts(3, 1), pts(3, 2), pts(3, 3)]);
        vecs(3, 2) = subs(Ast, [sr, st, sf] , [pts(3, 1), pts(3, 2), pts(3, 3)]);
        vecs(3, 3) = subs(Asf, [sr, st, sf] , [pts(3, 1), pts(3, 2), pts(3, 3)]);
        
        SRtrans = sphTOrec * AsS; % functions of sr, st, sf, rx, ry, rz
        fprintf('\nA_rx : %s \n', SRtrans(1, 1));
        vecs(1, 1) = subs(SRtrans(1, 1), [ sr, st, sf, rx, ry, rz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_ry : %s \n', SRtrans(2, 1));
        vecs(1, 2) = subs(SRtrans(2, 1), [ sr, st, sf, rx, ry, rz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_rz : %s \n', SRtrans(3, 1));
        vecs(1, 3) = subs(SRtrans(3, 1), [ sr, st, sf, rx, ry, rz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        
        SCtrans = sphTOcyn * AsS; % functions of sr, st, sf, cr, cf, cz
        fprintf('\nA_cr : %s \n', SCtrans(1, 1));
        vecs(2, 1) = subs(SCtrans(1, 1), [ sr, st, sf, cr, cf, cz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_cf : %s \n', SCtrans(2, 1));
        vecs(2, 2) = subs(SCtrans(2, 1), [ sr, st, sf, cr, cf, cz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] );
        fprintf('\nA_cz : %s \n', SCtrans(3, 1));
        vecs(2, 3) = subs(SCtrans(3, 1), [ sr, st, sf, cr, cf, cz ],...
            [ pts(3, 1), pts(3, 2), pts(3, 3), pts(1, 1), pts(1, 2), pts(1, 3) ] ); 
    end
    
    fprintf('\n    *** the vectors by component @ point *** \n');
    fprintf('rec :  %.4f  ,  %.4f  ,  %.4f  \n',...
        vecs(1, 1), vecs(1, 2), vecs(1, 3) );
    fprintf('cyn :  %.4f  ,  %.4f  ,  %.4f  \n',...
        vecs(2, 1), vecs(2, 2), vecs(2, 3) );
    fprintf('cph :  %.4f  ,  %.4f  ,  %.4f  \n',...
        vecs(3, 1), vecs(3, 2), vecs(3, 3) );

    fprintf('\n    *** all the magnitudes should be the same *** \n');
    fprintf('the magnitude of rec is: %.4f\n',...
        sqrt( vecs(1, 1)^2 + vecs(1, 2)^2 + vecs(1, 3)^2 ));
    fprintf('the magnitude of cyn is: %.4f\n',...
        sqrt( vecs(2, 1)^2 + vecs(2, 2)^2 + vecs(2, 3)^2 ));
    fprintf('the magnitude of sph is: %.4f\n',...
        sqrt( vecs(3, 1)^2 + vecs(3, 2)^2 + vecs(3, 3)^2 ));
end




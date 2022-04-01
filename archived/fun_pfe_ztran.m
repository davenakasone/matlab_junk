function zer_pol  = fun_pfe_ztran(fun_in)
%{
    make output whatever you need
        for zeros and poles:
        zer_pol  = fun_pfe_ztran(fun_in) ... zer_pol = [double(transpose(rotz_num)), double(transpose(rotz_den)) ];  

    you have to come back and change this for repeated poles

    sub out any finitie sums   !!! poly descend !!!
    make sure that fun_in has the proper form:
        
        num(zeros) :  (1 - (1/(Z * z1))) * (1 - (1/(Z * z2))) * ...  maybe * const
        den(poles) :  (1 - (1/(Z * p1))) * (1 - (1/(Z * p2))) * ...

        fun_in = num/den    in terms of "Z"
%}
global Z; % the running variable
fun_og = expand(fun_in); % make sure it is a descending polynomial(num и den)
tran_fun = 0;

    [num,den]=numden(fun_og);
    bcof = double(coeffs(num,Z,'all')); % numerator ...const absorbed
    acof = double(coeffs(den,Z,'all')); % denominator
    
    [ro, po, ko] = residuez(bcof, acof)
    len_ro = length(ro);
    len_po = length(po);
    len_ko = length(ko); 
    if len_ro > 0
        for ct = 1:len_ro
            temp = ro(ct,1)/(1-(1/Z)*po(ct,1));
            tran_fun = tran_fun + temp;
        end
    end
    
    if len_ko > 0
        for ct = 1:len_ko
            temp = Z^(1-ct) * ko(ct);
            tran_fun = tran_fun + temp;
        end
    end
    test1 = simplify(tran_fun-fun_in);

    [bi, ai] = residuez(ro, po, ko);
    len_bi = length(bi);
    len_ai = length(ai);
    temp_num = 0;
    temp_den = 0;
    for ct = 1:len_bi
        temp_num = temp_num + bi(ct)*Z^(len_bi-ct);
    end
    for ct = 1:len_ai
        temp_den = temp_den + ai(ct)*Z^(len_ai-ct);
    end
    temp_fun = temp_num / temp_den;
    test2 = simplify(temp_fun-fun_in);

    rotz_num = roots(bcof);
    rotz_den = roots(acof);
    t3n = 1;
    t3d = 1;
    if length(rotz_num) >= 1
        for ct = 1: length(rotz_num)
            t3n = t3n * ( Z - rotz_num(ct) );
        end
    end
    if length(rotz_den) >= 1
        for ct = 1:length(rotz_den)
            t3d = t3d * ( Z - rotz_den(ct) );
        end
    end
    t3 = t3n/t3d;
    test3 = simplify(t3-fun_in);
    
    zer_pol = [double(transpose(rotz_num)), double(transpose(rotz_den)) ];
    
    %
    fprintf('the original:\n');
    pretty(fun_in);
    fprintf('the PFE:\n');
    pretty(tran_fun);
    fprintf('zeros / poles form:\n');    % maybe off by a constant
    pretty(t3);
    fprintf('result1 is %s      ...0 means good\n', test1);
    fprintf('result2 is %s      ...0 means good\n', test2); 
    fprintf('result3 is %s      ...0 means good\n\n', test3);
    if length(rotz_num) < 1
        fprintf('there are no zeros\n');
    else
        fprintf('zeros, z = ');
        for ct = 1:length(rotz_num)
            fprintf('%d     ', rotz_num(ct,1));
        end
        fprintf('\n');
    end
    if length(rotz_den) < 1
        fprintf('there are no poles\n');
    else
        fprintf('poles, z = ');
        for ct = 1:length(rotz_den)
            fprintf('%d    ', rotz_den(ct,1));
        end
        fprintf('\n');
    end
    %}
end


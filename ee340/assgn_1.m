%{
    1  : problem 1
    2  : problem 2
    3  : problem 3
    4  : problem 3
    5  : problem 3 again
    6  : problem 3 with change    vpasolve()
%}
clc;
select = 6;

%-------------------------------------------------------------------------------------
if select == 1
    v_1 = 100 * exp(1j*deg2rad(0));
    v_2 = 100 * exp(1j*deg2rad(-120));
    v_3 = 100 * exp(1j*deg2rad(120));
    z_L = 0.2 + 1j;
    z_a = z_L + 1j*14 + 84;
    z_b = z_L + 1j*80 +70;
    z_c = z_L - (1j*18) + 45;
    
    
    % solve by kvl
    A = [z_a + z_b, -z_b     , 0;
         -z_b     , z_b + z_c, -z_c;
         0        , -z_c     , z_c + z_L];
    Y = [v_1 - v_2; v_2 - v_3; v_3];
    I = A \ Y;
    f_mdri("Ia", I(1,1), 1);
    f_mdri("Ib", I(2,1), 1);
    f_mdri("Ic", I(3,1), 1);
    fprintf("\n");
    i_1 = I(1,1);
    i_2 = I(2,1) - I(1,1);
    i_3 = I(3,1) - I(2,1);
    i_n = I(3,1);
    v_n = i_n * z_L;
    f_mdri("Vn", v_n, 1);
    f_mdri("In", i_n, 1);
    f_mdri("i1", i_1, 1);
    f_mdri("i2", i_2, 1);
    f_mdri("i3", i_3, 1);
    
    %check, by KCL
    fprintf("\n");
    syms vn;
    eqn1 = (vn - v_1) / z_a;
    eqn2 = (vn - v_2) / z_b;
    eqn3 = (vn - v_3) / z_c;
    eqn4 = (vn - 0) / z_L;
    kcl = eqn1 + eqn2 + eqn3 + eqn4;
    v_n = double(solve(kcl == 0, vn));
    i_n = v_n / z_L;
    i_1 = (v_1 - v_n) / z_a;
    i_2 = (v_2 - v_n) / z_b;
    i_3 = (v_3 - v_n) / z_c;
    f_mdri("Vn", v_n, 1);
    f_mdri("In", i_n, 1);
    f_mdri("i1", i_1, 1);
    f_mdri("i2", i_2, 1);
    f_mdri("i3", i_3, 1);
    %iii = i_1+i_2+i_3
    
    % power
    fprintf("\n");
    s_1 = v_1 * conj(i_1);
    s_2 = v_2 * conj(i_2);
    s_3 = v_3 * conj(i_3);
    s_n = v_n * conj(i_n);
    s_tot = s_1 + s_2 + s_3 + s_n;
    f_mdri("pow_1", s_1, 1);
    f_mdri("pow_2", s_2, 1);
    f_mdri("pow_3", s_3, 1);
    f_mdri("pow_n", s_n, 1);
    f_mdri("pow_tot", s_tot, 1);
    
    i2r_1 = abs(i_1)^2 * real(z_a);
    i2r_2 = abs(i_2)^2 * real(z_b);
    i2r_3 = abs(i_3)^2 * real(z_c);
    i2r_n = abs(i_n)^2 * real(z_L);
    check = i2r_1 + i2r_2 + i2r_3 + i2r_n
    checkk = real(s_tot) - check
end

%-------------------------------------------------------------------------------------
if select == 2
    v_1 = 120 * exp(1j*deg2rad(0));
    v_2 = 120 * exp(1j*deg2rad(-120));
    v_3 = 120 * exp(1j*deg2rad(120));
    R = 20;
    z_L = 1j*R*sqrt(3);
    z_C = -z_L;
    f_mdri("R", R, 1);
    f_mdri("z_L", z_L, 1);
    f_mdri("z_C", z_C, 1);
    fprintf("\n");
    
    v12 = v_1 - v_2;
    i12 = v12/R;
    v23 = v_2 - v_3;
    i23 = v23 / z_C;
    v31 = v_3 - v_1;
    i31 = v31 / z_L;
    f_mdri("v12", v12, 1);
    f_mdri("i12", i12, 1);
    f_mdri("v23", v23, 1);
    f_mdri("i23", i23, 1);
    f_mdri("v31", v31, 1);
    f_mdri("i31", i31, 1);
    fprintf("\n");
    
    I1 = i12 - i31;
    I2 = i23 - i12;
    I3 = i31 - i23;
    f_mdri("I2", I2, 1);
    f_mdri("I1", I1, 1);
    f_mdri("I3", I3, 1);
    fprintf("\n");
    
    p12 = v12 * conj(i12);
    p23 = v23 * conj(i23);
    p31 = v31 * conj(i31);
    f_mdri("p12", p12, 1);
    f_mdri("p23", p23, 1);
    f_mdri("p31", p31, 1);
    fprintf("\n");
    
    p1 = v_1 * conj(I1);
    p2 = v_2 * conj(I2);
    p3 = v_3 * conj(I3);
    f_mdri("p1", p1, 1);
    f_mdri("p2", p2, 1);
    f_mdri("p3", p3, 1);
    fprintf("\n");
    
    fprintf("\n\tso Q = 0, P = 720 W per phase, 2160 W total???\n");
end

%-------------------------------------------------------------------------------------
if select == 3
    z_line = 0.1 + 1j*0.2;
    S_mag  = 2300;
    pf = 0.8;
    tht = acos(pf);
    P = S_mag * pf;
    Q = S_mag * sin(tht);
    S = P + 1j*Q;
    
    syms V delta;
    v_s = 120 * V * exp(1j * delta);
    v_load = V * exp(1j * delta);
    I = (v_s - v_load) / z_line;
    
    eqn_L = z_line * conj(S);
    f_mdri("lhs", eqn_L, 1);
    eqn_R = 120*V - V^2 - eqn_L;
    vv = double(solve(eqn_R == 0, V));
    try_v1 = vv(1);
    try_v2 = vv(2);
    
    eqn_v1 = ( conj(S) / (try_v1 * exp(1j * delta) ) ) * z_line + try_v1 * exp(1j * delta);
    eqn_v2 = ( conj(S) / (try_v2 * exp(1j * delta) ) ) * z_line + try_v2 * exp(1j * delta);
    dd1 = double(solve( eqn_v1 == 120 * exp(1j * delta), delta))
    dd2 = double(solve( eqn_v2 == 120 * exp(1j * delta), delta))
    
    % group solution pairs and test all 4 combinations:
    sols_v = [try_v1, try_v2]
    sols_d = [dd1(1), dd2(1), 0];
    fprintf("\n\tpossible V solutions:  %0.3f + %0.3f  ,  %0.3f  %0.3f\n",...
        real(sols_v(1)), imag(sols_v(1)), real(sols_v(2)), imag(sols_v(2)));
    fprintf("\n\tpossible delta solutions, in radians:  %0.3f  ,  %0.3f  ,  %0.3f\n\n",...
        sols_d(1), sols_d(2), sols_d(3));
    for ii = 1:1:2
        for jj = 1:1:3
            temp = 120*exp(1j*sols_d(jj)) - (conj(S)/(sols_v(ii)*exp(1j*sols_d(jj))))*z_line - sols_v(ii)*exp(1j*sols_d(jj));
            fprintf("\t\tKCL difference:  real()=  %f  ,  imag()=  %f\n",...
                real(temp), imag(temp));
            temp = (120*exp(1j*sols_d(jj))-sols_v(ii)*exp(1j*sols_d(jj)))/z_line;
            chk_S = S - (sols_v(ii)*exp(1j*sols_d(jj))) * conj(temp);
            fprintf("\t\t\tdiff of S:  real()=  %f  ,  imag()=  %f\n\n",...
                real(chk_S), imag(chk_S));
        end
    end
    %chk_11 = 120*exp(1j*sols_d(1)) - (conj(S)/(sols_v(1)*exp(1j*sols_d(1))))*z_line - sols_v(1)*exp(1j*sols_d(1))
    %chk_12 = 120*exp(1j*sols_d(2)) - (conj(S)/(sols_v(1)*exp(1j*sols_d(2))))*z_line - sols_v(1)*exp(1j*sols_d(2))
    %chk_21 = 120*exp(1j*sols_d(1)) - (conj(S)/(sols_v(2)*exp(1j*sols_d(1))))*z_line - sols_v(2)*exp(1j*sols_d(1))
    %chk_12 = 120*exp(1j*sols_d(2)) - (conj(S)/(sols_v(2)*exp(1j*sols_d(2))))*z_line - sols_v(2)*exp(1j*sols_d(2))
end

%-------------------------------------------------------------------------------------
if select == 4
    z_line = 0.1 + 1j*0.2;
    S_mag  = 2300;
    pf = 0.8;
    tht = acos(pf);
    P = S_mag * pf;
    Q = S_mag * sin(tht);
    S = P + 1j*Q;
    
    syms V delta;
    v_s = 120 * V * exp(1j * delta);
    v_load = V * exp(1j * delta);
    
    rhs = simplify((v_s - v_load) * conj(v_load));
    lhs = z_line * conj(S);
end

%-------------------------------------------------------------------------------------
if select == 5
    z_line = 0.1 + 1j*0.2;
    S_mag = 2300;
    pf = 0.8;
    tht = acos(pf);
    P = S_mag * pf;
    Q = S_mag * sin(tht);
    S = P + 1j*Q;
        chk = sqrt(P^2 + Q^2);
        chkk = tht - angle(S);
    
    rhs = z_line * conj(S); % 460 + j230
    syms V;
    syms delta;
    eqn_real = 120 * V * cos(delta * V - delta) - real(rhs) - V^2;
    eqn_imag = 120 * V * sin(delta * V - delta) - imag(rhs);
    [sols_v, sols_d] = vpasolve([eqn_real == 0, eqn_imag == 0],[V, delta]);
    fprintf("\n\tV = %f  volts\n", sols_v);
    fprintf("\n\tdelta = %f  rad\n", sols_d);
    
    % check the work
    VV = sols_v;
    dd = sols_d;
    v_s = 120 * exp(1j*dd*VV);
    v_load = VV * exp(1j*dd);
    I = (v_s - v_load) / z_line;
    SS = v_load * conj(I);
    check = double(SS - S); % it is about 0, these values produce same power
end

%-------------------------------------------------------------------------------------
if select == 6
    z_line = 0.1 + 1j*0.2;
    S_mag = 2300;         % 2300 VA
    pf = 0.8;
    tht = acos(pf);
    P = S_mag * pf;       % 1840 W
    Q = S_mag * sin(tht); % 1380 var
    S = P + 1j*Q;
        chk = sqrt(P^2 + Q^2);
        chkk = tht - angle(S);
    
    rhs = z_line * conj(S); % 460 + j230
    syms delta;
    syms V;
    eq_r = 120 * V * cos(delta) - V^2 - real(rhs);
    eq_i = 120 * V * sin(delta) + imag(rhs);
    [dd, vv] = vpasolve([eq_r==0, eq_i==0], [delta, V]);
    
    % verify
    fprintf("delta = %0.3f  rad  ,  V = %f  volts\n", dd, vv);
    v_load = vv * exp(1j*dd);
    I = (v_s - v_load) / z_line;
    KVL = v_s - I * z_line - v_load;
    SS = v_s*conj(I) - (I*z_line)*conj(I);
    pff = cos(angle(SS));
    fprintf("\tKVL differenece:  real= %f  ,  image= %f\n", real(KVL), imag(KVL));
    fprintf("\t\tS: pf= %0.2f  ,  abs()= %f  ,  real()= %f  ,  imag()= %f\n\n",...
        pff, abs(SS), real(SS), imag(SS)); 
    fprintf("\nnormalizing delta to 2 pi interval, delta = %f  radians\n", dd+(36*2*pi));
end
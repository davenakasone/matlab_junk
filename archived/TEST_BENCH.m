%{
                            TEST BENCH

    const = cls_CONST('helper');   % create your obj to get physical constants
    ee = cls_EE330_helper();       % create your obj for EE330 calculations

    fun_conT(yourFun, inpR, pnts)....continous time
        % needs function in terms pt, specify range, any points, get a 2D graph

    fun_conT2(funA, inpA, funB, inpB, pnts)
        % 2 functions, one figure
    
    fun_conTabg(yourFun, alpha, beta, gama, inpR, pnts)
        % same as fun_conT , but shows adjusted function

    fun_conTx(funs, range)
        all on seperate graphs

    fun_CT_zero_pol(zers,pols, roc, side)
        give it the roots and it will make a pol/zero for the Laplace in s-domain
    
    rootz = fun_demovFun()
        % manual approach to get roots form specific formula

    [myRoots] = fun_demovNum(numIn,num, den)    
        % for getting all the roots of a complex radical

    fun_disN(yourFun, inpR, pnts)
        % simple discrete function global in
    fun_disN2(funA, inpA, funB, inpB, pnts)
    fun_disN3(funA, inpA, funB, inpB, funC, inpC, pnts)
    fun_disNstr8(inpt, dpnt, title_nm)   only numbers
    fun_disNx(funs, range)   all on same screen
    fun_disNxxx(funs, range)  seperate screens

    fun_DT_zer_pol(zers,pols, roc, side)
        pole-zero plot for Z-trans

    fun_graph2Dcomplex(pointsIn)
        % just for plotting complex points 2D , put in as [ 1 x n ] pt = z = x + jy,  prin 1st

    fun_graph_cyn(axis,rad, cent, hgt, pts)
        % graph a cylinder and any rect points

    fun_graph_norm(funIn, pts)
        % give it f(rx, ry) and it gives funIn as surface with normal vectors on it

    fun_graph_path3(pts, Vf)
        % take a path of points, and a vector field to see a line integral

    fun_graph_spaceC(params, range)      can even do complex
        % give it [ x(t), y(t), z(t) ] and a range for t    use global param pt;

    fun_graph_sph(rad, cent, pts)
        % graph a sphere and any rect points

    fun_graphEllipsoid(cent,semAx)
        % just needs [x, y, z] center and length of axis...everything else gets ellipsoid

    fun_graphSurCon(funIn, select)        select: 1, 2, 3, 4, 5, 6              
        % big 4: real, imag, abs, phase + a contour plot(real vs imag)    

    fun_graphSurCon2(funIn1, funIn2, select)   select: 1, 2, 3, 4, 5, 6    funz1 surf, funz2 mesh
        % same as single funtion version, but 2 functions

    fun_graphSurConX(funs, select)
        % multiple figures

    fun_graphVF(field, pts, levZ)
        % supposed to take a vector field and quiver it      and see it 2D with levlZ

    result = fun_harminic2test(funIn)
        % 2D laplace check    result 0 on fail, 1 on pass 

    tran_fun = fun_pfe_ztran(fun_in, checks)
        give it Z-transform in proper form, makes a PFE for easy inverse Z trans

    result = fun_rCRtest(funIn)         result = 1 pass, result = 0 fail        
        % Cauchy-Riemen 2 part test 

    fun_tetraView(funIn, centReal, centImag)    % center on point, shows input and output plane

    publisher('filename') ...pdf, html, ect

    include global variables (used as symbols) as specified below
    don't force a global to take a fixed value (mutate is bad)
    don't use a naked 'j' as a variable     1j equiv sqrt(-1) or try to index with it
    don't use reserves like pi, ans, and other confusing shit
    don't overload "s" и "t" while using Laplace and Fourier built-in functions

TODO list:

    - validate feed() to only accept points [0, pi]   ...input control on all the functions of eeHelper
    - proof all the spots atan() could cause problems in eeHelper...use atan2() if possible
    - control /by 0 for graphing funcitons
    - fun_graphsurCon2() could be improved to (fun1, fun2, sel1, sel2) ...lots of ifs 6x6 opts
    - need an esp safety for all the graphs ... no /by 0
    - figure('Name', 'xxxx', ...) can help get all the fun_graphXX()'s to talk to eachother global ax?
    - fundamental transform derivations
    - camelCase or use_under_scores ???

    - use more LateX
    - get those cell blocks working for different size vectors {}
%}
%color = uisetcolor([1, 1, 0], 'Selecf Color');    % [.9, .9, .9] is a nice gray
%clear all;
clc;
close all;
clearvars;
%sympref('MatrixWithSquareBrackets', 1);
sympref('PolynomialDisplayStyle', 'descend');   % usually 'descend' is best...  or ascend
old_val = sympref('HeavisideAtOrigin', 1);
format shortE; % default, short, long, shortE, longE, shortG, longG, +, hex, rational 
format compact; % [compact,loose]

                % try not to mutate the globals too much

% Universal constants                
global null; null = double.empty();  % nulled array of size 1x1
global intv; intv = eps*(1e14);  % a interval of  0.0222   to get around /by 0, impulse, ect
global j; j = sqrt(-1); % just overload it here so no problems occur later

% Greeks
global alp; syms alp; assume(alp, 'real'); % real phase, real z
global bet; syms bet; assume(bet, 'real'); % imag phase, imag z
global chie; syms chie; assume(chie, 'real'); % polarization
global chim; syms chim; assume(chim, 'real'); % magnetization
global ep; syms ep; assume(ep, 'real'); % as in ep = ep0 * epr    or D = ep * E
global ep0; syms ep0; assume(ep0, 'real'); % sub with const.ep0   permittivity of free space
global eta; syms eta; % intrinisic impedance for coupling Maxwell's system of equations
global gam; syms gam; % phase coeff
global lam; syms lam; assume(lam, 'real'); % wave length, eigen val, 
global mu; syms mu; assume(mu, 'real'); % as in mu = mu0 * mur    or B = mu * H , friction coeff
global mu0; syms mu0; assume(mu0, 'real'); %  sub with const.mu0   permiability of free space
global sig; syms sig; assume(sig, 'real'); % resistivity , mean, 
global omg; syms omg; assume(omg, 'real'); % omega , angular freq ... 2 pi f
global phi; syms phi; assume(phi, 'real'); % phase constant
global tau; syms tau; assume(tau, 'real'); % convolution dummy, torque, trans coeff
global tht; syms tht; assume(tht, 'real'); % any angle theta
global zta; syms zta; assume(zta, 'real'); % any zeta

% Single Vars (specific)
global amp; syms amp; assume(amp, 'real'); % amplitude of any function
global freq; syms freq; assume(freq, 'real'); % in Hz  f = 1 / T  = 2 pi / omg
global phs; syms phs; assume(phs, 'real'); % phasor term, (omg*t - beta*rz), pi/2, ect
global c1; syms c1; % integration constant 1
global c2; syms c2; % integration constant 2
global c3; syms c3; % integration constant 3
global c4; syms c4; % integration constant 4
global pL; syms pL; assume(pL, 'real'); assume(pL >= 0); % rho_line charge denisty in C / m
global pS; syms pS; assume(pS, 'real'); assume(pS >= 0); % rho_surface charge denisty in C / m^2
global pV; syms pV; assume(pV, 'real'); assume(pV >= 0); % rho_volume charge denisty in C / m^3

% Single Vars (general)
global k; syms k; assume(k, {'real', 'integer'}); % index k
global n; syms n; assume(n, {'real', 'integer'}); % index n, iztrans
global n0; syms n0; assume(n0, {'real', 'integer'}); % arbitrary n0, usually offset
global N; syms N; assume(N, {'real', 'integer'});  % bound of series, DT period
global r; syms r; assume(r, 'real'); %  z = r exp(j tht)
global s; syms s; % j omg, laplace
global t; syms t; assume(t,'real'); % time , ilaplace
global t0; syms t0; assume(t0, 'real'); % arbitrary t0 offset
global T; syms T; assume(T, 'real');  % CT period, integration limit
global U; syms U; assume(U, 'real');   % U( X, Y)  ...f(Zxy) = U(X,Y) + j V(X,Y)  --> real(f)     
global V; syms V; assume(V, 'real');   % V( X, Y)  ...f(Zxy) = U(X,Y) + j V(X,Y)  --> imag(f)
global X; syms X; assume(X, 'real'); % real Z
global Y; syms Y; assume(Y, 'real'); % imag Z
global Z; syms Z;   % do not compound, z-trans
global Zxy; syms Zxy; Zxy = X + 1j*Y; % compound at will

% Coords       phi --> fi
global rx; global ry; global rz; % "rectangular x"   "rectangular y"  "rectangular z"       
global cr; global cf; global cz; % "cylindrical rho" "cylindrical fi" "cylindrical z"       
global sr; global st; global sf; % "spherical radius" "spherical theta" "spherical fi"   
syms rx; assume(rx, 'real'); syms ry; assume(ry,'real'); syms rz; assume(rz, 'real');
syms cr; assume(cr, 'real'); assume(cr >= 0);
syms cf; assume(cf, 'real'); assume(cf >= 0);
syms cz; assume(cz, 'real');
syms sr; assume(sr, 'real'); assume(sr >= 0);
syms st; assume(st, 'real'); assume(st >= 0);
syms sf; assume(sf, 'real'); assume(sf >= 0);

global arx; syms arx; assume( arx, 'real'); % x_hat (rec)
global ary; syms ary; assume( ary, 'real'); % y_hat (rec)
global arz; syms arz; assume( arz, 'real'); % z_hat (rec)
global acr; syms acr; assume( acr, 'real'); % r_hat (cyn)
global acf; syms acf; assume( acf, 'real'); % phi_hat (cyn)
global acz; syms acz; assume( acz, 'real'); % z_hat (cyn)
global asr; syms asr; assume( asr, 'real'); % r_hat (sph)
global ast; syms ast; assume( ast, 'real'); % theta_hat (sph)
global asf; syms asf; assume( asf, 'real'); % phi_hat (sph)

global Arx; syms Arx; assume( Arx, 'real'); % x_hat comp (rec)
global Ary; syms Ary; assume( Ary, 'real'); % y_hat comp (rec)
global Arz; syms Arz; assume( Arz, 'real'); % z_hat comp (rec)
global Acr; syms Acr; assume( Acr, 'real'); assume(Acr >= 0); % r_hat comp (cyn)
global Acf; syms Acf; assume( Acf, 'real'); assume(Acf >= 0); % phi_hat comp (cyn)
global Acz; syms Acz; assume( Acz, 'real'); % z_hat comp (cyn)
global Asr; syms Asr; assume( Asr, 'real'); assume(Asr >= 0); % r_hat comp (sph)
global Ast; syms Ast; assume( Ast, 'real'); assume(Ast >= 0);% theta_hat comp (sph)
global Asf; syms Asf; assume( Asf, 'real'); assume(Asf >= 0);% phi_hat comp (sph)


ee = cls_EE330_helper();
const = cls_CONST();
mat = cls_math('Themistocles');
%publisher('ee330_hw4.m');


                            select = 0;  % CHANGE CHANGE CHANGE

%fprintf("\nselection #%d\n\n", select);
%------------------------------------------------------------------------------------------ #0
if select == 0 
    inpR = [-1, 0, 1, 2, 3];
    r
end


%------------------------------------------------------------------------------------------ #1111
if select == 1111
    
    
end


%------------------------------------------------------------------------------------------ #2222
if select == 2222
    display('function test');
end


%------------------------------------------------------------------------------------------ #3333
if select == 3333
    clc;
    display('full circle conversion');
    pt_rec = [-1, 1, -1];                            % CHANGE
    ee.feed(pt_rec, 'R');
    %ee.print();
    pt_cyn = ee.pts(2,:);
    pt_sph = ee.pts(3,:);
    
    
    vf_rec = [ rx, ry, rz];                          % CHANGE
    fprintf('vf_rec:   [ %s , %s , %s ]  --> cyn\n', vf_rec); % RC -> CR
    vf_cyn = ee.transVecRC(vf_rec);
    fprintf('vf_cyn:   [ %s , %s , %s ]   --> rec\n', vf_cyn);  
    vf_rec1 = ee.transVecCR(vf_cyn);
    fprintf('vf_rec1:  [ %s , %s , %s ]  --> MATCH\n', simplify(vf_rec1));
    
    fprintf('\nvf_rec:   [ %s , %s , %s ]  --> sph\n', vf_rec); % RS -> SR
    vf_sph = ee.transVecRS(vf_rec);
    fprintf('vf_sph:   [ %s , %s , %s ]    --> rec\n', vf_sph);
    vf_rec2 = ee.transVecSR(vf_sph);
    fprintf('vf_rec2:  [ %s , %s , %s ]  --> MATCH\n', simplify(vf_rec2));
    
    toss = input('any key:', 's');
    clc;
    
    vf_cyn = [sin(cf)*cr^2, cos(cf)*cz, cz^2];   % CHANGE
    fprintf('\nvf_cyn:   [ %s , %s , %s ] --> rec\n', vf_cyn);    % CR -> RC
    vf_rec = ee.transVecCR(vf_cyn);
    fprintf('vf_rec:   [ %s , %s , %s ] --> cyn\n', simplify(vf_rec));
    vf_cyn1 = ee.transVecRC(vf_rec);
    fprintf('vf_cyn1:  [ %s , %s , %s ] --> MATCH\n', simplify(vf_cyn1));
    
    fprintf('\nvf_cyn:   [ %s , %s , %s ] --> sph\n', vf_cyn);    % CS -> SC
    vf_sph = ee.transVecCS(vf_cyn);
    fprintf('vf_sph:   [ %s , %s , %s ] --> cyn\n', simplify(vf_sph));
    vf_cyn2 = ee.transVecSC(vf_sph);
    fprintf('vf_cyn2:  [ %s , %s , %s ] --> MATCH\n', simplify(vf_cyn2));
    
    toss = input('any key:', 's');
    clc;
    
    vf_sph = [ sr^2, cos(st)*sin(sf)*sr, cos(sf)*sr];  % CHANGE
    fprintf('\nvf_sph:   [ %s , %s , %s ] --> rec\n', vf_sph);    % SR -> RS
    vf_rec = ee.transVecSR(vf_sph);                             
    fprintf('vf_rec:   [ %s , %s , %s ] --> sph\n', simplify(vf_rec));
    vf_sph1 = ee.transVecRS(vf_rec);
    fprintf('vf_sph1:  [ %s , %s , %s ] --> MATCH\n', simplify(vf_sph1));
    
    fprintf('\nvf_sph:   [ %s , %s , %s ] --> cyn\n', vf_sph);    % SC -> CS
    vf_cyn = ee.transVecSC(vf_sph);
    fprintf('vf_cyn:   [ %s , %s , %s ] --> sph\n', simplify(vf_cyn));
    vf_sph2 = ee.transVecCS(vf_cyn);
    fprintf('vf_sph2:  [ %s , %s , %s ] --> MATCH\n', simplify(vf_sph2));
end


%------------------------------------------------------------------------------------------ #4444
if select == 4444
    cm.funz = Z^2;
    pretty(cm.funz);
    fprintf('%s\n',cm.name);
    
    tst = cm.test_static_ex(1);
    fprintf('%d\n',tst);
    cm.test_static();
    
    tst = cm.test_self(-1);
    fprintf('%d\n', tst);
    tst = cm.test_self_ex(2);
    fprintf('%d\n', tst);
end

fprintf("\n\n\t\t~ ~ ~ PROGRAM COMPLETE ~ ~ ~\n\n");

%------------------------------------------------------------------------------------------
%******************************************************************************************
%------------------------------------------------------------------------------------------
%global fig3r; fig3r = figure();    that is how you plot to same graph...just pass fig handle
%burner(fig3r);

%{
                             PROPOGANDA
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.^ .*  ./   ...   \  '  .'   @       %%                   \delta[n]  LaTex
%.4f  %s  ...
==   &&  &  |  ||    ~=

arry = ClassName.empty(); if need empty array

abs()
angle()
bode()
break;
cart2sph()
cat(dim, A, B)
charpoly()
coeffs(fun,var,'All') get coeffs of poly...
collect()
compass()
conj()
conv()
cross(A, B)
cumsum([1, 2, 3]) --> [1, 3, 6]
deg2rad()
delete()  clear out a syms or handle
diag()
diff(fun, var, order)  
 dirac()   unit impulse
divergence(X, Y, Z, Fx, Fy, Fz)...not good cyn и sph
drawnow
double()
dot(A, B)
dsolve()
eig()
end ...array(2:2:end) skips every other
eps
expand()
eye()  ...identity matrix
factor()
factorial()
find( arr<0)
filter()
fix()  round to 0
fliplr()
flipup()
fminbnd()
fminsearch()
freqs()
for, while, continue, break, if, else, else if
fourier()
fzero()  can get roots of almsot anything, not just polys
getframe
global
gradient(scalar fun, [vars])
heaviside()... 1/2 at var = 0
imag()
inf
ifourier()
IgnoreAnalyticConstraints
ilaplace(F,var,tranVar) usually (F, s, t)
impulse()
Inf
input(prompt, 's')
int(fun, var, low, high)   symbolic --> subs
intmax('uint64')
integral()  all numeric
inv()   inverse of square matrix...use A\B for sys...
isempty()
isinf()
isnan()
isolate(eqn,xpr)
iztrans()
kron()
laplace(f,var,tranVar) usually (f,t,s)
length()  largest array dimension
limit()
linespace(min, max, points)
log()  base e
log2() base 2
log10()  base 10
max(arr, [], 'all')
meshgrid( dim1, dim2, dim3)  or just 2D
min(arr, [], 'all')  'linear'
mod()   b = mod(a,m)   b = mod(23,5) = 3
NaN
nchoosek(n,k)   = factorial(n) / ( factorial(n-k) * factorial(k) )
no-op
norm(A)
num2str()
[N,D] = numden(A)   get numerator and denominator
numel()
rad2deg()
real()  
[r,p,k] = residue(b,a)  good for partial fractions
[b,a] = residue(r,p,k)
rewrite()
round()
roots()
sum()
sys
switch case
optimset()
partfrac()  .... pf = partfrac(f,z,'FactorMode', 'full')   ...get that pfe done right
pascal()
pause(secs)
perms()
piecewise()
pol2cart()
poly()  inverse of roots()
polyder()
polyfit()
polyint()
polyval()
potential()
pretty(yourEqn)
rand()
randi()
randperm()
reallog()  no neg or complex
realmax()
residue()
residuez()
rng()
simplify()
simplifyFraction()   ...pfe partial fraction expansion
[rows, cols,...] = size(array);
solve()      solves almost anything
sph2cart()
sprintf()
subs(fun, old, new)
sympref()
syms vs sym() ...
symsum()
taylor()
tf()
[A,B,C,D] = tf2ss(b,a)
tf2zp()
transpose()
vectorPotential()
vid = videoWriter('movName', 'MPEG-4')  , open(vid) , writeVideo(vid, flipBook)  , close(vid) see spaceC
ztrans()
zp2tf()


quiver, mesh, meshc surf, fsurf, cylinder, sphere, fplot3(space curve), contour   ...camlight 
line, scatter, fplot()...fcontour, fimplicit, fmesh, fsurf,...
zplane(), bode(), fplot, area(), patch()

silly trig and complex exp
f = exp(afaf) cosfafaf...abortion
f = rewrite(f, 'sincos')
Ref = simplify(collect(real(f)))        % Real part
Imf = simplify(collect(imag(f)))        % Imaginary part


%------------------------------------------------------------------------------------------ 
>> help ops
  Operators and special characters.
 
  Arithmetic operators.
    plus       - Plus                               +    
    uplus      - Unary plus                         +    
    minus      - Minus                              -    
    uminus     - Unary minus                        -    
    mtimes     - Matrix multiply                    *    
    times      - Array multiply                    .*    
    mpower     - Matrix power                       ^    
    power      - Array power                       .^    
    mldivide   - Backslash or left matrix divide    \    
    mrdivide   - Slash or right matrix divide       /    
    ldivide    - Left array divide                 .\    
    rdivide    - Right array divide                ./    
    idivide    - Integer division with rounding option.
    kron       - Kronecker tensor product
    pagemtimes - Page-wise matrix multiply
    pagetranspose  - Page-wise transpose
    pagectranspose - Page-wise complex conjugate transpose
 
  Relational operators.
    eq         - Equal                             ==     
    ne         - Not equal                         ~=     
    lt         - Less than                          <      
    gt         - Greater than                       >      
    le         - Less than or equal                <=     
    ge         - Greater than or equal             >=     
 
  Logical operators.
    relop      - Short-circuit logical AND         &&     
    relop      - Short-circuit logical OR          ||     
    and        - Element-wise logical AND           &      
    or         - Element-wise logical OR            |      
    not        - Logical NOT                        ~      
    punct      - Ignore function argument or output ~
    xor        - Logical EXCLUSIVE OR
    any        - True if any element of vector is nonzero
    all        - True if all elements of vector are nonzero
 
  Special characters. 
    colon      - Colon                              : 
    paren      - Parentheses and subscripting      ( )              
    paren      - Brackets                          [ ]     
    paren      - Braces and subscripting           { }          
    punct      - Function handle creation           @
    punct      - Decimal point                      .      
    punct      - Structure field access             .      
    punct      - Parent directory                   ..     
    punct      - Continuation                       ...    
    punct      - Separator                          ,      
    punct      - Semicolon                          ;      
    punct      - Comment                            %      
    punct      - Invoke operating system command    !            
    punct      - Assignment                         =
    punct      - Quote                              '      
    punct      - Double quote                       "    
    transpose  - Transpose                         .'
    ctranspose - Complex conjugate transpose        ' 
    horzcat    - Horizontal concatenation          [,]     
    vertcat    - Vertical concatenation            [;]     
    subsasgn   - Subscripted assignment          ( ),{ },.   
    subsref    - Subscripted reference           ( ),{ },.   
    numArgumentsFromSubscript - Number of arguments for indexing methods
    subsindex  - Subscript index
    metaclass  - Metaclass for MATLAB class         ?
 
  Bitwise operators.
    bitand     - Bit-wise AND.
    bitcmp     - Complement bits.
    bitor      - Bit-wise OR.
    bitxor     - Bit-wise XOR.
    bitset     - Set bit.
    bitget     - Get bit.
    bitshift   - Bit-wise shift.
 
  Set operators.
    union      - Set union.
    unique     - Set unique.
    intersect  - Set intersection.
    setdiff    - Set difference.
    setxor     - Set exclusive-or.
    ismember   - True for set member.
%}


%{
    % poly long division     divisor|----  quo, rem      dividend (num) / divisor (den)   
                                      %dividend              f = quo + rem/divisor
    dividend = 2*s^3 + 7*s^2 + 4*s +9; 
    coef_divd = double(coeffs(dividend, s, 'all'));  % [ 2, 7, 4, 9]
    divisor = s^2 + 1;
    coef_divs = double(coeffs(divisor, s, 'all')); % [ 1, 0, 1 ]
    [quo, rem] = deconv(coef_divd, coef_divs); % quo= 2x + 7 [2,7] : rem = 2x + 2 [0,0,2,2]
    check = expand(((2*s + 7) * divisor) + (2*s + 2)); %      quo*divisor + rem = dividend
%}


%{
    % Z-trans PFE   using  ex 10.9, p758   make your terms (1-1/zp)  ...  (1-1/zz)

    z1 = 3-(5/(6*Z));
    p1 = 1-(1/(4*Z));
    p2 = 1-(1/(3*Z));
    xz = z1 / (p1*p2);
    xz = expand(xz);
    [num,den]=numden(xz);
    bcof = double(coeffs(num,Z,'all'));
    acof = double(coeffs(den,Z,'all'));
    [ro, po, ko] = residuez(bcof, acof); % leave this open to plug in, or make function
    t1 = ro(1,1)/(1-(1/Z)*po(1,1));
    t2 = ro(2,1)/(1-(1/Z)*po(2,1));
    temp = t1+t2;
    check1 = simplify(temp-xz); % should be 0
    [bi, ai] = residuez(ro, po, ko); % leave open to compose
    xnum = bi(1)*Z^2 + bi(2)*Z^1 + bi(3)*Z^0;
    xden = ai(1)*Z^2 + ai(2)*Z^1 + ai(3)*Z^0;
    x_comp = xnum/xden;
    check2 = simplify(x_comp-xz); % should be 0   ...you brought it full circle

    tst = fun_pfe_ztran(xz, 1); ..good user-defined function
%}

%{
    poly division

    u = [ 2 , 7 , 4 , 9 ]; % 2 x^3 + 7 x^2 + 4 x + 9
    v = [ 1, 0, 1 ];       % x^2 + 1
    [q,r] = deconv(u,v);    v|u  quo, remainder    check  quo + rem/v = u
%}

%{
    PFE
    
    b=expand((s+1)*(2*s+11));         % 2s^2 + 13s +11
    b = double(coeffs(b,s,'all'));
    a=expand((6+5*s+s^2)*(s+4));      % s^3 + 9s^2 +26s + 24
    a = double(coeffs(a,s,'all'));
    [R,P,K] = residue(b,a);  %  or [b,a] = residue(r,p,k)   
%}

%{
    Cauchy Reimann Analytic Test:

    f = cos(Z);
    fz = subs(f, Z, Zxy);
    fun_rCRtest(fz);
%}

%{
    general limit

    f = (2*Z-5*sin(Z^2)) / (3*Z);
    res = limit(f, Z, Inf);
%}

%{
    Integration of an analytic function in simply connected domain  ... 0 if no holes and loop

    start = 1;
    stop = 2;
    f = sin(Z);
    F = int(f, Z);
    integral = subs(F, Z, stop) - subs(F, Z, start);
%}

%{
    Integration of non-analytic by parameterization
    
    fz = 1/Z;

    path1 = cos(pt) + 1j*sin(pt); % change at will based on contour
    start1 = 0;
    stop1 = 2*sym(pi);
    d_path1 = diff(path1, pt, 1);
    intg1 = subs(fz, Z, path1) * d_path1;
    intgg1 = int(intg1, pt);
    ul1 = subs(intgg1, pt, stop1);
    ll1 = subs(intgg1, pt, start1)
    intggg1 = ul1 - ll1;

    path2 = cos(pt) + 1j*sin(pt); % change at will based on contour
    start2 = 0;
    stop2 = 2*sym(pi);
    d_path2 = diff(path2, pt, 1);
    intg2 = subs(fz, Z, path2) * d_path2;
    intgg2 = int(intg2, pt);
    ul2 = subs(intgg2, pt, stop2);
    ll2 = subs(intgg2, pt, start2)
    intggg2 = ul2 - ll2;
%}

%{
    Taylor:
    
    z0 = 3;  % center 
    f = log(Z);  % function to Taylor
    a0 = (1/factorial(0))*subs(f, Z, z0); % first coeff a0
    
    f1 = diff(f, Z, 1);
    a1 = (1/factorial(1))*subs(f1, Z, z0);
    
    f2 = diff(f, Z, 2);
    a2 = (1/factorial(2))*subs(f2, Z, z0);
    
    f3 = diff(f, Z, 3);
    a3 = (1/factorial(3))*subs(f3, Z, z0);
    
    f4 = diff(f, Z, 4);
    a4 = (1/factorial(4))*subs(f4, Z, z0);
    
    f5 = diff(f, Z, 5);
    a5 = (1/factorial(5))*subs(f5, Z, z0);

    f6 = diff(f, Z, 6);
    a6 = (1/factorial(6))*subs(f6, Z, z0);
    
    f7 = diff(f, Z, 7);
    a7 = (1/factorial(7))*subs(f7, Z, z0);
%}

%{
    Binomial expansion (radical or large negative exponent downstairs
                    (1+z)^k   form
    format rational;
    k = -1/2;
    z0 = 1j;
    term = -1/4;  % less Z^n , if not 1
    const = 1/2;

    a0 = ( const * term^0 * 1 );
    a0 = a0 / factorial(0);

    a1 = ( const * term^1 * (k) );
    a1 = a1 / factorial(1);

    a2 = ( const * term^2 * (k) * (k-1) );
    a2 = a2 / factorial(2);

    a3 = ( const * term^3 * (k) * (k-1) * (k-2) );
    a3 = a3 / factorial(3); 

    a4 = ( const * term^4 * (k) * (k-1) * (k-2) * (k-3) );
    a4 = a4 / factorial(4); 

    a5 = ( const * term^5 * (k) * (k-1) * (k-2) * (k-3) * (k-4) );
    a5 = a5 / factorial(5);
    
    a6 = ( const * term^6 * (k) * (k-1) * (k-2) * (k-3) * (k-4) * (k-5) );
    a6 = a6 / factorial(6);
%}

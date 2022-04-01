clc;

syms kk;
syms nn;
syms pp;

pdf_bin(kk, nn, pp) = (factorial(nn)/(factorial(kk) * factorial(nn-kk)))*(pp.^kk)*(1-pp).^(nn-kk);

success = 1;     % k
trials = 3;      % n
prob = 0.5;     % p

fprintf("\n\tp =  %0.4f    ;    n = %3d    ;    k =  %3d\n", prob, trials, success);
quickk = pdf_bin(success, trials, prob);
fprintf("\n\t\tP(x == %2d) :  %0.4f\n", success, quickk);

summer = 0;
temp = 0;
fprintf("\n\t  X   P(X = k)    P(X <= k)\n");
for ii = 0:1:trials
    temp = pdf_bin(ii, trials, prob);
    summer = summer + temp;
    fprintf("\t %2d    %0.4f       %0.4f\n", ii, temp, summer);
end


%
%lcm(6,8)

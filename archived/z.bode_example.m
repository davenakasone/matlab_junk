%{

H(s) = 200s / (s^2 + 12s + 20)

%}

numerator = [200, 0];
denominator = [1, 12, 20];
bode(numerator, denominator);

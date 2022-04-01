%homo
syms t;
syms y1h(t);
syms y2h(t);
syms y3h(t);
A = [11, 6, 18; 9, 8, 18; -9, -6, -16]; % change me

Yh = [y1h(t); y2h(t); y3h(t)];
odesH = diff(Yh) == A*Yh;

[sy1h(t), sy2h(t), sy3h(t)] = dsolve(odesH);
sy1h(t) = simplify(sy1h(t));
sy2h(t) = simplify(sy2h(t));
sy3h(t) = simplify(sy3h(t));

homo = [ sy1h(t); sy2h(t); sy3h(t) ];
display(homo);



%particular
syms y1p(t);
syms y2p(t);
syms y3p(t);
B = [0; 3; 4]; % change me

Yp = [y1p(t); y2p(t); y3p(t)];
odesP = diff(Yp) == A*Yp + B; % change me

[sy1p(t), sy2p(t), sy3p(t)] = dsolve(odesP);
sy1p(t) = simplify(sy1p(t));
sy2p(t) = simplify(sy2p(t));
sy3p(t) = simplify(sy3p(t));

parituclar = simplify([ sy1p(t); sy2p(t); sy3p(t) ]- homo); % change me a little....after - , keep 0 or - homo
display(parituclar);




%ivp
syms y1(t);
syms y2(t);
syms y3(t);

Y = [y1(t); y2(t); y3(t)];
odes = diff(Y) == A*Y + B; % change me
cond1 = y1(0) == -2; % change me
cond2 = y2(0) == 4; % change me
cond3 = y3(0) == 3; % change me
cons = [cond1; cond2; cond3];

[sy1(t), sy2(t), sy3(t)] = dsolve(odes, cons);
sy1(t) = simplify(sy1(t));
sy2(t) = simplify(sy2(t));
sy3(t) = simplify(sy3(t));

combined = [ sy1(t); sy2(t); sy3(t) ];
display(combined);


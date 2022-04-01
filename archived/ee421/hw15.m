%{
    1  : p11.11a, Vol
    2  : p11.11a, Vsp
    3  : p11.11b, Voh
    4  : p11.11b, Vsp

%}
clc;
selector = 4;

%---------------------------------------------------------------------------------------
if (selector == 1)
    syms V;
    lam = (30e3) * (120e-6) * (10/1);
    eqn = lam * (4.2 * V - (1/2) * V^2) + V;
    Vol = solve( eqn == 5, V);
    fprintf("\t\tlow output voltage:  %0.0f mV\n\n", Vol*1000);
    double(subs(eqn, V, Vol(1)))
end


%---------------------------------------------------------------------------------------
if (selector == 2)
    syms V;
    lam = (30e3) * (120e-6) * (10/1);
    eqn = lam * ((V - 0.8) * V - (1/2) * V.^2) + V;
    Vsp = solve( eqn == 5, V);
    fprintf("\t\tswitching point voltage:  %0.3f V\n\n", Vsp);
    double(subs(eqn, V, Vsp(2)))
end


%---------------------------------------------------------------------------------------
if (selector == 3)
    syms V;
    lam = (50e3) * (40e-6) * (10/1);
    eqn = lam * ((4.1) * (5 - V) - (1/2) * (5 - V).^2 ) - V;
    Voh = solve(eqn == 0, V);
    fprintf("\t\toutput high voltage:  %0.3f V\n\n", Voh);
    double(subs(eqn, V, Voh(1)))
end


%---------------------------------------------------------------------------------------
if (selector == 4)
    syms V;
    lam = (50e3) * (40e-6) * (10/1);
    eqn = lam * ((5 - V - 0.9) * (5 - V) - (1/2) * (5 - V).^2 ) - V;
    Vsp = solve(eqn == 0, V);
    fprintf("\t\tswitching point voltage:  %0.3f V\n\n", Vsp);
    double(subs(eqn, V, Vsp(1)))
end


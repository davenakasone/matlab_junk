function lhrh = fun_lft(zs, ws)
    global Z;
    syms z;
    syms w;
    z1 = zs(1); z2 = zs(2); z3 = zs(3);
    w1 = ws(1); w2 = ws(2); w3 = ws(3);
    lhs = ((w-w1)/(w-w3)) * ((w2-w3)/(w2-w1));
    rhs = ((z-z1)/(z-z3)) * ((z2-z3)/(z2-z1));
    %lhrh = [lhs, rhs];
    
    gam = ((Z-z1)/(Z-z3)) * ((z2-z3)/(z2-z1)) * ((w2-w1)/(w2-w3));
    solz = (w1-w3*gam)/(1-gam);
    fprintf('w(z)=\n');
    pretty(simplify(solz));
    lhrh = solz;
    
    alp = ((w-w1)/(w-w3)) * ((w2-w3)/(w2-w1)) * ((z2-z1)/(z2-z3));
    solw = (z3*alp)/(1-alp);
    fprintf('z(w)=\n');
    pretty(simplify(solw));
end


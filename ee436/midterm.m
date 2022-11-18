%{
    exam1

    1: demo TE_10 with known solutions
    2: p2e, find alpha_c
    3: p3, integrating for fields
    4: p3, demonstrate boundary condition
%}
clc;
close;
clear all;
select = 3;


if select == 2
    syms A_; %assume(A_, 'real'); assumeAlso(A_ > 0);
    syms a; %assume(a, 'real'); assumeAlso(a > 0);
    syms b; %assume(b, 'real'); assumeAlso(b > 0);
    syms beta;
    syms kc;
    syms mu; %assume(mu, 'real'); assumeAlso(mu > 0);
    syms Rs; %assume(Rs, 'real'); assumeAlso(Rs > 0);
    syms w; %assume(w, 'real'); assumeAlso(w > 0);
    syms x;
    syms y;
    m = 2;
    n = 0;
    E_x = (1j * w * mu * n * sym(pi) / (b * kc^2)) * A_ * cos(m * sym(pi) * x / a) * sin(n * sym(pi) * y / b);
    E_y = (-1j * w * mu * m * sym(pi) / (a * kc^2)) * A_ * sin(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    E_z = 0;
    H_x = (1j * beta * m * sym(pi) / (a * kc^2)) * A_ * sin(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    H_y = (1j * beta * n * sym(pi) / (b * kc^2)) * A_ * cos(m * sym(pi) * x / a) * sin(n * sym(pi) * y / b);
    H_z = A_ * cos(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    fprintf("\nthe E field:\n\n");
    E = [E_x, E_y, E_z];
    pretty(E);
    fprintf("\n the H field:\n\n");
    H = [H_x, H_y, H_z];
    pretty(H);

    % power down guide, no conduction
    fprintf("\npoynting vector:\n\n");
    p = cross(E, conj(H));
    pretty(p);
    igrand = dot(p, [0,0,1]);
    fprintf("\nintegrand is ready:\n\n");
    pretty(igrand);
    intg_y = int(igrand, y, [0, b]);
    fprintf("\nafter integrating y[0:b] ->\n\n");
    pretty(intg_y);
    intg_x = int(intg_y, x, [0, a]);
    integrated = (1/2) * intg_x;
    fprintf("\nafter integrating x[0:a], P_%d%d ->\n\n", m, n);
    pretty(integrated);
    fprintf("\nsubstituting kc = sqrt([m pi / a]^2 + [n pi / b]^2), P_%d%d :\n\n", m, n);
    P_ = subs(integrated, kc, sqrt((m * sym(pi) / a)^2 + (n * sym(pi) / b)^2));
    pretty(P_);

    % setting to the TE_20 form + dimensions given
    fprintf("\nsubstituting b==a/3, P_%d%d :\n\n", m, n);
    P20 = subs(P_, b, a/3);
    pretty(P20);

    % power loss in a cross section
    Js_y = cross([1, 0, 0], [0, 0, 1]) .* subs(H_z, x, 0); % x == 0  wall
    Js_z = cross([0, 1, 0], [1, 0, 0]) .* subs(H_x, y, 0); % y == 0 wall
    Js_x = cross([0, 1, 0], [0, 0, 1]) .* subs(H_z, y, 0); % y == 0 wall
    Js = Js_x + Js_y + Js_z;
    fprintf("\nJs:\n\n");
    pretty(Js);
    intg_x0 = dot(Js_y, Js_y);
    int_x0 = int(intg_x0, y, [0, b]); % Js_y, x == 0 wall
    intg_y0a = dot(Js_x, Js_x);
    int_y0a = int(intg_y0a,  x, [0, a]); % Js_x, y == 0 wall
    intg_y0b = dot(Js_z, Js_z);
    int_y0b = int(intg_y0b,  x, [0, a]); % Js_z, y == 0 wall
    fprintf("\nPl:\n\n");
    Pl = simplify(subs((int_x0 + int_y0a + int_y0b) * Rs, kc, sqrt((m * sym(pi) / a)^2 + (n * sym(pi) / b)^2)));
    pretty(Pl);

    % setting to the TE_20 form + dimensions given
    fprintf("\nsubstituting b==a/3, P_l :\n\n", m, n);
    P_l = subs(Pl, b, a/3);
    pretty(P_l);

    % now find alpha_c
    alpha_c = P_l / (2 * P20);
    fprintf("\nalpha_c:\n\n");
    pretty(alpha_c);
end


%%%%~~~~


if select == 3
    syms A_; 
    syms a;
    syms beta;
    syms kc;
    syms m;
    syms mu; 
    syms n;
    syms w; 
    syms x;
    syms y;
    
    H_z = A_ * (cos(m * sym(pi) * x / a) * cos(n * sym(pi) * y / a) +...
                cos(n * sym(pi) * x / a) * cos(m * sym(pi) * y / a));
    H_x = (-1j * beta / kc^2) * diff(H_z, x);
    H_y = (-1j * beta / kc^2) * diff(H_z, y);
    H = [H_x, H_y, H_z];

    E_x = (-1j * w * mu / kc^2) * diff(H_z, y);
    E_y = (1j * w * mu / kc^2) * diff(H_z, x);
    E_z = 0;
    E = [E_x, E_y, E_z];

    fprintf("\nHx:\n\n");
    pretty(H_x);
    fprintf("\nHy:\n\n");
    pretty(H_y);
    fprintf("\nHz:\n\n");
    pretty(H_z);

    fprintf("\nEx:\n\n");
    pretty(E_x);
    fprintf("\nEy:\n\n");
    pretty(E_y);
    fprintf("\nEz== %d\n\n", E_z); 
end


%%%%~~~~


if select == 4
    syms A_; 
    syms a;
    syms beta;
    syms kc;
    syms m;
    syms mu; 
    syms n;
    syms w; 
    syms x;
    syms y;
    
    H_z = A_ * (cos(m * sym(pi) * x / a) * cos(n * sym(pi) * y / a) +...
                cos(n * sym(pi) * x / a) * cos(m * sym(pi) * y / a));
    H_x = (-1j * beta / kc^2) * diff(H_z, x);
    H_y = (-1j * beta / kc^2) * diff(H_z, y);
    H = [H_x, H_y, H_z];

    E_x = (-1j * w * mu / kc^2) * diff(H_z, y);
    E_y = (1j * w * mu / kc^2) * diff(H_z, x);
    E_z = 0;
    E = [E_x, E_y, E_z];

    Hz_dx = diff(H_z, x);
    Hz_dy = diff(H_z, y);
    Hz_dx_sub = subs(Hz_dx, y, x);
    Hz_dy_sub = subs(Hz_dy, y, x);
    eqn = Hz_dx_sub - Hz_dy_sub;

    fprintf("\nd/dx{Hz}:\n\n");
    pretty(Hz_dx);
    fprintf("\nd/dy{Hz}:\n\n");
    pretty(Hz_dy);

    fprintf("\nd/dx{Hz} y = x :\n\n");
    pretty(Hz_dx_sub);
    fprintf("\nd/dy{Hz} y = x :\n\n");
    pretty(Hz_dy_sub);

    fprintf("\nd/dy{Hz} -  d/dy{Hz}   with:: y = a-x -->\n\n");
    pretty(simplify(eqn));
end


%%%%~~~~


%%%%~~~~END>  midterm.m


if select == 99
    fprintf("\n\n\t\tdone\n\n");
end
































if select == 1
    syms A_; %assume(A_, 'real'); assumeAlso(A_ > 0);
    syms a; %assume(a, 'real'); assumeAlso(a > 0);
    syms b; %assume(b, 'real'); assumeAlso(b > 0);
    syms beta;
    syms kc;
    syms mu; %assume(mu, 'real'); assumeAlso(mu > 0);
    syms Rs; %assume(Rs, 'real'); assumeAlso(Rs > 0);
    syms w; %assume(w, 'real'); assumeAlso(w > 0);
    syms x;
    syms y;
    m = 1;
    n = 0;
    E_x = (1j * w * mu * n * sym(pi) / (b * kc^2)) * A_ * cos(m * sym(pi) * x / a) * sin(n * sym(pi) * y / b);
    E_y = (-1j * w * mu * m * sym(pi) / (a * kc^2)) * A_ * sin(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    E_z = 0;
    H_x = (1j * beta * m * sym(pi) / (a * kc^2)) * A_ * sin(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    H_y = (1j * beta * n * sym(pi) / (b * kc^2)) * A_ * cos(m * sym(pi) * x / a) * sin(n * sym(pi) * y / b);
    H_z = A_ * cos(m * sym(pi) * x / a) * cos(n * sym(pi) * y / b);
    fprintf("\nthe E field:\n\n");
    E = [E_x, E_y, E_z];
    pretty(E);
    fprintf("\n the H field:\n\n");
    H = [H_x, H_y, H_z];
    pretty(H);
    
    % power down guide, no conduction
    fprintf("\npoynting vector:\n\n");
    p = cross(E, conj(H));
    pretty(p);
    igrand = dot(p, [0,0,1]);
    fprintf("\nintegrand is ready:\n\n");
    pretty(igrand);
    intg_y = int(igrand, y, [0, b]);
    fprintf("\nafter integrating y[0:b] ->\n\n");
    pretty(intg_y);
    intg_x = int(intg_y, x, [0, a]);
    integrated = (1/2) * intg_x;
    fprintf("\nafter integrating x[0:a], P_%d%d ->\n\n", m, n);
    pretty(integrated);
    fprintf("\nsubstituting kc = sqrt([m pi / a]^2 + [n pi / b]^2), P_%d%d :\n\n", m, n);
    P_ = subs(integrated, kc, sqrt((m * sym(pi) / a)^2 + (n * sym(pi) / b)^2));
    pretty(P_);

    % power loss in a cross section
    Js_y = cross([1, 0, 0], [0, 0, 1]) .* subs(H_z, x, 0); % x == 0  wall
    Js_z = cross([0, 1, 0], [1, 0, 0]) .* subs(H_x, y, 0); % y == 0 wall
    Js_x = cross([0, 1, 0], [0, 0, 1]) .* subs(H_z, y, 0); % y == 0 wall
    Js = Js_x + Js_y + Js_z;
    fprintf("\nJs:\n\n");
    pretty(Js);
    intg_x0 = dot(Js_y, Js_y);
    int_x0 = int(intg_x0, y, [0, b]); % Js_y, x == 0 wall
    intg_y0a = dot(Js_x, Js_x);
    int_y0a = int(intg_y0a,  x, [0, a]); % Js_x, y == 0 wall
    intg_y0b = dot(Js_z, Js_z);
    int_y0b = int(intg_y0b,  x, [0, a]); % Js_z, y == 0 wall
    fprintf("\nPl:\n\n");
    Pl = simplify(subs((int_x0 + int_y0a + int_y0b) * Rs, kc, sqrt((m * sym(pi) / a)^2 + (n * sym(pi) / b)^2)));
    pretty(Pl);

    % now find alpha_c
    alpha_c = Pl / (2 * P_);
    fprintf("\nalpha_c:\n\n");
    pretty(alpha_c);
end
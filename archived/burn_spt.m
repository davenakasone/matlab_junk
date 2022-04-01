clc;
syms sig;
syms omg;
pts = 10;

sigF = linspace(-10,10,pts);
omgF = sigF;

[SIG, OMG]= meshgrid(sigF,omgF);

fun = sig^2 + omg^2;
F = zeros(pts,pts);

for row = 1:pts
    for col = 1:pts
        F(row,col) = SIG(row,col)^2 + OMG(row,col)^2;
    end
end

figure();
colormap jet;
colorbar;
meshc(SIG,OMG,F);

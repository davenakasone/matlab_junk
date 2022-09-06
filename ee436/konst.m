classdef konst
    properties (Constant)
        c = 299.8e6;             % speed of light, m/s
        mu0 = 4*sym(pi)/10^7;    % permeability of free space, H/m
        ep0 = 8.854e-10;         % permittivity of free space, F/m
        eta0 = 377;              % wave impedance in free space, ohms
    end
end
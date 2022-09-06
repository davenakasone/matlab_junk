function val_perm = f_perm_mu(mu_r)
%{
    given relative permeability, multiplies by permeability of free space
    if you want mu_0, just make argument "1"
    permeability measures ability of material to allow magnetic lines of force to pass
%}
    mu_0 = 4 * pi / 10^7; % in H / m
    val_perm = mu_r * mu_0;
end
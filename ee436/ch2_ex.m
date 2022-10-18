%{
    chapter 1, EM

    1   :  try smith chart

%}
clc;
clear;
select = 1;


if select == 1
   Z_load = 40 + 1j*70;
   Z_0 = 100;

   z_load = Z_load / Z_0;
   p_cplx("z_load", z_load);

   G = (z_load - 1) / (z_load + 1);
   p_cplx("Gamma", G);
end


%%%%~~~~
if select == 1
   
end


%%%%%%%%~~~~~~~~END>  ch2_ex.m

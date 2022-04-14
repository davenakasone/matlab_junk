function z_out = f_z_motor(R1, X1, Xm, R2, X2, s_s)
    temp1 = R1 + 1j*X1;
    temp2 = (R2/s_s) + 1j*X2;
    temp3 = f_para(1j*Xm, temp2);
    z_out = temp1 + temp3;
end
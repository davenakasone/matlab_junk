function f_rec2pol(rect_in)
    in_abs = abs(rect_in);
    in_ang = angle(rect_in);
    fprintf("abs():  %0.3f  ,  angle():  %0.3f\n",...
        in_abs, rad2deg(in_ang));
end
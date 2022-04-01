%{

    rectangular/cartesian points are structured   (recx, recy, recz)         GREEN o
    rectangular/cartesian vectors are structured < recx, recy, recz >
    recx = x , recy = y , recz = z

    cylindrical/circular points are structured (cynr, cynf, cynz)            BLUE  x
    cylindrical/circular points are structured < cynr, cynf, cynz >
    cynr = rho  ,  cynf = fi  ,  cynz = z

    spherical points are structured (sphr, spht, sphf)                       RED   *
    spherical points are structured < sphr, spht, sphf > 
    sphr = r  ,  spht = theta  ,  sphf = fi

    designed to take input as a { vector or point }   treated same... just a tip trace
    specify rect (rectangular)  , cyn (cylindrical) , or sph (spherical)
    output will be a 3x3 with all the representations
    
         rec   cyn   sph          diagonals are redundant, you will get them based on input
    rec   x     1     1
    cyn   1     x     1
    sph   1     1     x 



    input points 
  
        comp1   comp2   comp3           this has a rectangular input and needs cyn и sph      
   rec    a1     a1      a1
   cyn    x      x       x
   sph    x      x       x 



    eventually, output will be
    
        comp1   comp2   comp3                 
   rec    a1     a1      a3       < a1, a2, a3 >      ==      < recx, recy, recz >
   cyn    b2     b2      b3       < b1, b2, b3 >      ==      < cynr, cynf, cynz >
   sph    c3     c3      c3       < c1, c2, c3 >      ==      < sphr, spht, sphf >

    use them however you want
        see check with cart2sphy....
%}
clc; 
clf;
close all;
clearvars;
select = 13; 
%color = uisetcolor([1, 1, 0], 'Select Color'); % .9, .9, .9 is nice

points = NaN(3,3);
input = [1, 2, 3];     % CHANGE HERE         *** only inputs of entire script ***
type = 'rec';          %  only use 'rect'  'cyn'  or 'sph'

if type == 'rec'
    points(1,1) = input(1);
    points(1,2) = input(2);
    points(1,3) = input(3);
end
if type == 'cyn'
    points(2,1) = input(1);
    points(2,2) = input(2);
    points(2,3) = input(3);
end
if type == 'sph'
    points(3,1) = input(1);
    points(3,2) = input(2);
    points(3,3) = input(3);
end





syms x;
x = 2;
display(points);

display('fuck you');
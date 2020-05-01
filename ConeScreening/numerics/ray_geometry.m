clear all
c = 343;
y = [0.1:0.1:1];
x0 = 1;
l_s = 0.8;

a  = sqrt(x0^2 + l_s^2+ y.^2);
alpha = 0.0456;

r = sqrt((1+0.8*cos(alpha)^2)+ (0.8*alpha+y).^2);

dlambda = 2*(r-a);

nu = c./dlambda









function p = ellipseperimeter(r1,r2)
%
% Use the Ramanujan formula to approximate the perimeter of an ellipse
%
% https://www.mathsisfun.com/geometry/ellipse-perimeter.html
% Approximation 2
%
p = pi * (3 * (r1 + r2) - sqrt(3*r1 + r2) * (r1 + 3*r2));
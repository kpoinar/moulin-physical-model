function p = ellipseperimeter(r1,r2)
%
% Use the Ramanujan formula to approximate the perimeter of an ellipse
%
% https://www.mathsisfun.com/geometry/ellipse-perimeter.html
% Approximation 2
%
% CT July 16, 2020, simplified the function so that 
% it only output the perimeter and not both the perimeter and the Area

% r1: the major axis
% r2: the minor axis

p = pi.* (3 .*(r1 + r2) - sqrt((3.* r1 + r2) .* (r1 +3 .* r2))); 

%a = pi .* r1 .* r2;
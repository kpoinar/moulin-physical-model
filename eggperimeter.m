function p = eggperimeter(r1,r2)
%
% Our moulin is an egg: half ellipse, half circle
%
p = 0.5*ellipseperimeter(r1,r2) + pi*r1;
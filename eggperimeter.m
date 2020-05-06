function p = eggperimeter(r1,r2)
%
% Our moulin is an egg: half ellipse, half circle
[ep, ~] = ellipseperimeter(r1, r2);
p = 0.5*ep + pi*r1;
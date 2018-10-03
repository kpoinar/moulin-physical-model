function u = conserveWaterMass(Mr,z,u0,z0)
%
% Given the water velocity u0 at position z0,
% and a known moulin radius at each position z,
% calculate the water flux at each height z.
%
% Find index of the given water velocity at z0
[~,jj] = min(abs(z-z0));
%
% Calcualte u at all other points in z
u = u0*Mr(jj)^2 ./ Mr.^2;
u(jj) = u0;   % that should already be correct

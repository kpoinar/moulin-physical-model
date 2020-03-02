function dV = waterVolumeFromTurbulence(r1, r2, dr1, dr2, z, wet, dt)
%
% Compare the old area and the new area at each level z
% Area of moulin:
% half circle + half ellipse
% A = 1/2 pi*r1^2 + 1/2 pi * r1 * r2
%   = pi/2 * r1 * (r1 + r2)
%
C = makeConstants;
%
Aold = pi/2 * r1 .* (r1 + r2);
Anew = pi/2 * (r1 + dr1) .* (r1 + r2 + dr1 + dr2);
%
% Area change:
dA = Anew - Aold;
%
% Sum over all areas that have water
dV = trapz(z(wet),dA(wet));
%
% That is ice volume.  Convert to water volume (it is lower):
dV = dV * C.rhoi / C.rhow;
%
% And divide by dt?
dV = dV / dt;


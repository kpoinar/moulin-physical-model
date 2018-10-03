function [dM,Vadd] = turbulence(Mr,u,L,dt,H,hw,z,C)
%
% Turbulent melting is a function of moulin radius (Mr), water speed (u), 
% and hydraulic gradient (Psi).
%
% Lauren derived this for a cylinder using Schoof 2010, who modeled a
% half-cylinder subglacial channel.
%
% Hydraulic gradient:
Psi = C.rhow*C.g* hw / L;
% Psi(z>hw) = 0;
% Psi = C.rhow*C.g*H/L;
%
% factor f:
f = 4/C.rhow*Mr.*Psi ./u.^2;



%f = f * 1/100;






%f = 1;%0.1;
% constants (from Schoof 2010):
c1 = 1/C.rhoi/C.Lf;
c4 = sqrt(2*pi^3/C.rhow./f);
%
% change (growth) in moulin radius from turbulent melting per timestep:
% dM = dt / (2*pi*C.Lf) * Mr .* Psi .* sqrt(2*pi*Mr.*Psi./C.rhow./f);
dM = dt/2/pi*c1*c4 .* Mr.^1.5 .* Psi.^1.5;
dM(z>hw) = 0;
%
% How much water volume is added?
Vadd = C.rhoi/C.rhow * 2*pi*trapz(z,Mr.*dM);
%
% How does this affect moulin radius?
%Mr = Mr + dM;
%

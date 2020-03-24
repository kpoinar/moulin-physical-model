% Catania and Neumann 2010, Equation 1

function [dF, Vfrz] = refreeze_simple(Mrminor,Mrmajor,Tfar,z,hw,t,dt,C)

%Constant = 2093; % J/kg/K %(CT)I don't know this constant. What should we call it?
%kappa = 37.2/31557600; % thermal diffusivity, m2/s : 1.18e-6 (CT)replaced
%this value with kappa in makeConstant, but value is slightly different.
% KP - This is the specific heat capacity of ice, using the value from Alley (2005), http://doi.org/10.3189/172756405781813483 .  We have our own value of ice specific heat capacity in makeConstants.m, as C.cp = 2115 J/kg/K.  We should use C.cp here.
dF = C.cp*(Tfar-T0)*2*sqrt(C.kappa*t/pi)/C.Lf;
dFprev = C.cp*(Tfar-T0)*2*sqrt(C.kappa*(t-dt)/pi)/C.Lf;
dF = dF + dFprev;
% No freezing where there's no water
dF(z>hw) = 0;
% ^ The value of dF is negative (closure) because Tfar < To
%
% Calculate new moulin radius, while making sure that dF didn't freeze 
% more than the moulin diameter
Mrminornew = max(Mrminor + dF,0);
Mrmajornew = max(Mrmajor + dF,0);
dF_major = Mrmajornew - Mrmajor; % implement ceiling
dF_minor = Mrminornew - Mrminor; % implement ceiling


% How much water loss occurred to refreezing?
Vfrz = C.rhoi/C.rhow * (trapz(z,pi*Mrminor.*dF_minor) + trapz(z,0.5*ellipseperimeter(Mrminor,Mrmajor).*dF_major));
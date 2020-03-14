% Catania and Neumann 2010, Equation 1

function [Mr, dF] = refreeze_simple(Mr,Tfar,z,hw,t,dt,C)

Constant = 2093; % J/kg/K %(CT)I don't know this constant. What should we call it?
%kappa = 37.2/31557600; % thermal diffusivity, m2/s : 1.18e-6 (CT)replaced
%this value with kappa in makeConstant, but value is slightly different.

dF = Constant*(Tfar-T0)*2*sqrt(C.kappa*t/pi)/C.Lf;
dFprev = Constant*(Tfar-T0)*2*sqrt(C.kappa*(t-dt)/pi)/C.Lf;
dF = dF + dFprev;
% No freezing where there's no water
dF(z>hw) = 0;
% ^ The value of dF is negative (closure) because Tfar < To
%
% Calculate new moulin radius, while making sure that dF didn't freeze 
% more than the moulin diameter
Mrnew = max(Mr + dF,0);
dF = Mrnew - Mr; %implement ceiling
Mr = Mrnew;
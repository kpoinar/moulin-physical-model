% Catania and Neumann 2010, Equation 1

function [Mr, dF] = refreeze_simple(Mr,Tfar,z,hw,t,dt)

C = 2093; % J/kg/K
Lf = 335000; % J/kg
kappa = 37.2/31557600; % thermal diffusivity, m2/s : 1.18e-6
dF = C*(Tfar-273.15)*2*sqrt(kappa*t/pi)/Lf;
dFprev = C*(Tfar-273.15)*2*sqrt(kappa*(t-dt)/pi)/Lf;
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
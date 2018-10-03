% Based on thermal crevasse model (but 1D, not 2D)
% 
%
    % Next improvement will be with varying water levels in the borehole,
    % which will vary the boundary condition.  Sometimes a Stefan problem 
    % (no longer with a continuous time of exposure), other times a simple
    % conduction problem.
    % Colin Meyer's course notes were helpful for understanding the Stefan 
    % problem:  http://people.seas.harvard.edu/~colinrmeyer/Solidification%20of%20Fluids%20Notes.pdf
    %
    % A problem here will be switching BACK to a Stefan solution after some
    % time away from it.  During the time that the water level in the
    % moulin is low, and the ice in question is exposed to cold air
    % temperatures, it will reset to a cool temperature near the interface,
    % thus losing some time with respect to the continuous Stefan solution.
    % How to treat?
    % Compare the temperature profile (Ttrue) in the ice (by solving heat 
    % equation at each time step, with Lf going in as a source term at the 
    % bdy?) to the hypothetical temperature profile (Thyp) if the ice were 
    % continously in contact with the water, and find the time t* for which 
    % Thyp(t*) = Ttrue(t), within some tolerance (0.01°C maybe) at some 
    % distance (1 m?) away from the boundary.  Then continue using the 
    % Stefan solution, but with t* instead of t. 
%
%
function [Mr, dF, Tnew, Vfrz] = refreeze(Mr,Tzx,z,hw,dFprev,nx,x,dx,dt)
%
% Do 1D or 2D here?
% Steepest temperature gradient across 50 meters in the vertical on the
% Luthi profile: ~0.5°C / 50 meters = 0.01°C/m
% Expected temperature gradient across the 80 meter domain in x: 
% Tfar (Luthi) reaches temperature minimum of -22°C
% Twater is 0°C, expect this to be ~50-100 m away (depending how long I
% choose to make the x domain)
% Thus, dTdx ~ 0.3°C/m at the coldest point in the Luthi bow
% dTdx will match dTdz ~ 0.01°C only within the bottom 2% of the z domain,
% where T(z) > 0.8°C.
% This area is small AND unimportant for refreezing, since the ice is
% nearly temperate close to the bed.
% So I can exclude vertical diffusion from the model, run it in 1D only.
%
% C = 2093; % J/kg/K
% Lf = 335000; % J/kg
% kappa = 37.2/31557600; % thermal diffusivity, m2/s : 1.18e-6
% dF = C*(Tfar-273.15)*2*sqrt(kappa*t/pi)/Lf;
% dFprev = C*(Tfar-273.15)*2*sqrt(kappa*(t-dt)/pi)/Lf;
% dF = dF + dFprev;
%
% Solve for ice temperature:
Tnew = NaN(size(Tzx));
% Holder for the temperature gradient near the moulin wall:
Txscale = NaN(size(Tzx,1),2);
% Constants
ki =    2.1; % J/mKs
rhoi = 910; % kg/m3
rhow = 1000; % kg/m3
cp = 2115;  % J/kgK
Lf = 3.35e5; % J/kg
C = ki /(rhoi * cp); 
% The diffusion lengthscale:
xscale = sqrt(C * dt);
for r = 1:size(Tzx,1)
    Tnew(r,:) = solveTmoulin(Tzx(r,:)',Tzx(r,1),Tzx(r,end),C,dFprev(r),nx,dx',dt);
    Txscale(r,:) = interp1(x,Tnew(r,:),[0 xscale]);
end
% 
% Translate the ice temperature gradient into a refreezing thickness for
% this timestep
dF = -abs(dt * ki / rhoi / Lf*diff(Txscale,[],2)/xscale);
%
% No freezing where there's no water
dF(z>hw) = 0;
% ^ The value of dF is negative (closure) because Tfar < To
%
% Calculate new moulin radius, while making sure that dF didn't freeze 
% more than the moulin diameter
Mrnew = max(Mr + dF,0);
dF = Mrnew - Mr;
Mr = Mrnew;
%
% How much refroze?
Vfrz = rhoi/rhow * trapz(z,2*pi*Mr.*dF);
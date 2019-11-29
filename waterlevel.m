function [hw, cumvol, M0, min_index] = waterlevel(Mr,z,V0)
% Solve for water level in the moulin based on the known volume of water we
% have to work with, and the moulin radius (size) at all levels z
%
% Cumvol = cumtrapz(z,pi*Mr.^2);
% j = find(Cumvol < V0,1,'last');
% hw = z(j);

% lca completed it differently
cumvol          = cumtrapz(z,pi*Mr.^2);
[M0, min_index] = min(abs(V0-cumvol));
extra_water     =  V0-cumvol(min_index);
hw              = z(min_index) + extra_water ./(pi .* (Mr(min_index).^2) );
hw(hw<0) = 0;

%
%
%
%
% Don't let it get too low
%hw = max(hw, 0.1*z(end));
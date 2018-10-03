function hw = waterlevel(Mr,z,V0)
% Solve for water level in the moulin based on the known volume of water we
% have to work with, and the moulin radius (size) at all levels z
%
Cumvol = cumtrapz(z,pi*Mr.^2);
j = find(Cumvol < V0,1,'last');
hw = z(j);
%
%
%
%
%
% Don't let it get too low
hw = max(hw, 0.1*z(end));
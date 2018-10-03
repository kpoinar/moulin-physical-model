function [Qout, Qc, Qd, Sc, hd] =subglacial(hw, H, dt, Sc, dx, dy, C, ffsub, k0)
%
% Moulin discharge (Qout) is dependent on the height of the water within the
% moulin (hw) and the capacity of the subglacial system to accept water.
% hw = moulin water level (m)
% H  = ice thickness (m)
% dt = time step (s)
% Sc = channel cross-sectional area
% dx = distance between terminus and moulin
% dy = width of the distributed system (very clearly user defined)
% C  = typical constants rhow rhoi g, etc
% ffsub = friction factor for the subglacial system - because this is a
% lumped element model, the subglacial ff will be constant at 650
% kg/m^(8/3)
% k0 = dimensionless transmissivity coefficient 10^(-6)
% C.nw = need to add the water viscosity to the constants -> 10^(-3) Pa s 
% C.ni = ice viscosity -  can be calculated (see Hoffman and Price 2014
% equation 17)
% This module is developed using the the less complicated formulation of
% the a channel-distributed system from Hoffman and Price and Hewitt.
% Because the needs of the moulin model are relatively simple (it just
% needs the volume of water accepted by the subglacial system, this
% formulation is similar to a lumped-element analysis.
%
%
%% Distributed sheet evolution 
% I think that there will need to be a pde solver


%% Channel evolution

% 1. calculate discharge through the channel (Qc) using the water level in
% the moulin
Qc = (1 / ffsub)^(1/2) .* Sc^(4/3) .* C.rhow .* C.g .* hw;

% 2. calculate the evolution of the subglacial system
% channel melting
Mc  = (1 / C.Lf) * (Qc * C.rhow * C.g * hw) * dx; % The hydraulic gradient 
%used here assumes that dx is the distance from the terminus to the moulin
%and the terminus water level is effectively 0.

% 3. Calculate mass balance of water in the system

% exchange between channel and distributed system, which makes the
% assumption that the far edge of the distributed system has a water level
% equivalent to the fraction of overburden
lambda = - (k0 * hd^3 / C.nw) * (C.rhow * C.g * hw -  C.rhow * C.g * 0.917*H) / dy; 

dSc_water =  dt * ((Mc / rhow) * dx + lambda * dx - Qc);

% 4. Evolve the channel area
%calculate the effective pressure of the channel 
Nc = 
dSc_ice  = Mc / C.rhoi - Sc


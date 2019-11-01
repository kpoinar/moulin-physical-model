function [dM, uw, Vadd] =turbulence_hutter(hw, Qout, Mr, z, dt)
%
%!!! note, currently there is no mechanism to deal with temperature, so
%include_ice_temperature should be set as false!!!!
% Moulin discharge (Qout) is dependent on the height of the water within the
% moulin (hw) and the capacity of the subglacial system to accept water.

%inputs from the moulin model
%hw   = hight of the water in the moulin (m)
%Qout = discharge of water into the subglacial system (m3/s)
%Mr   = current moulin radius (m)
%Ti   = ice temperature in *C (from Luthi import file)
%dt   = timestep  (s)
%z    = elevation of each model node 0 = bottom of moulin; H= z(top of
%moulin)
%Tw    = water temperature in the moulin 
C = makeConstants;

%exports
%dM  = change in moulin radius due to melting (m) 
%uw  = water velocity at each node (m/s)
%Vadd = volume of water added due to melting for each grid cell (m3)

% internal variables
manrough = 0.03;
fr = 0.1;

include_ice_temperature = false; %True means that the melt is partially dependent 
                                %on the ice and water temperature (e.g. Clarke 2002) 
                                %In this case the water temperature must evolve
                                %via several equations. False
                                %means that melting is only dependent on
                                %turbulence (e.g. Gully et al., 2014)

if include_ice_temperature
    Ti = importTz('Luthi', z);
else
    Ti = NaN;
end


%initialize other variables


%% Initiate calculations for turbulent melting 

% logical matrix to determine in what nodes water is present (head given)
% or not present ==0;
waterpresent = hw - z;
waterpresent(waterpresent<0) =0;


uw = Qout ./Mr;  %calculate water velocity within the moulin column
uw(waterpresent ==0) =0; % if there is no water in a given cell,


% calculate moulin cross-sectional area from radius
S = pi .* Mr .^2; 

% calculate the effective pressure 
Pi   = C.rhoi .* C.g .* flipud(z); % ice pressure
Pw   = C.rhow .*C.g .* waterpresent;

Mp   = 2 .* pi .* Mr; % wetted/melting perimeter
Rh   = Mr ./2; %hydraulic radius 
Dh   = (4*(pi .* Mr.^2)) ./ Mp; %hydrualic diameter

Re   = 4 .* C.rhow .* abs(uw) .* Rh  ./ C.mu; %reynolds number
Pr   = C.mu .* C.cp ./ C.kw;
Nu   = 0.023 .* Re .^(4/5) .* Pr .^(2/5);

fR   = 8 .* C.g .* (manrough.^2) ./ (Rh.^(1/3)); %Darcy weisbach friction factor for varying conduit geometry
tau0 = (1/8) .* fR .* C.rhoi .* uw .* abs(uw); % wall stress exerted by turbulent flow 


if include_ice_temperature
    
melt   = (Mp .* C.kw .* Nu .* (Tw - Ti) ./ (4 .* C.Lf .* Rh)) ./( dz);
   % now the units are the same for both the melts...

dTwdz = (diff(Tw) ./ diff(z));
dTwdz  = [0, dTw];
mdot = melt .* ( dz);

dTw =  -uw .* dTwdz  + 1 ./ (C.rhow .* C.cp .* S) ...
       .* (Mp .* tau0 .* uw - mdot .* (C.Lf + C.cp .* ...
       (Tw - Ti) - ((uw.^2)./2)));
   
else
dz         = nanmean(diff(z)); %find the length of the nodes to use as the length scale
head_loss  = ((uw.^2) .* fr .* dz) ./(2.* Dh .* C.g);
melt       = (Qout .* C.rhow .* C.g .* (head_loss ./ dz))./ ...
             (2 .* pi .* C.Lf .* Mr); % these units are m per second of wall melt
         
end

dM  = melt .* dt; %change in radius over the given time step
Vadd = (C.rhoi ./ C.rhow) .*(pi .* ((Mr+dM) - Mr)); %volume of meltwater for each node



end







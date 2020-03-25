function [dM, uw, Vadd, hL] =turbulence_hutter(hw, Qout, Mr, z, dt)

%
% Moulin discharge (Qout) is dependent on the height of the water within the
% moulin (hw) and the capacity of the subglacial system to accept water.

%inputs from the moulin model
%hw   = height of the water in the moulin (m)
%Qout = discharge of water into the subglacial system (m3/s)
%Mr   = current moulin radius (m)
%Ti   = ice temperature in *C (from Luthi import file)
%Tw   = water temperature in the moulin
%z    = elevation of each model node 0 = bottom of moulin; H= z(top of
%moulin)
%dt   = timestep  (s)
C = makeConstants;

%exports
%dM  = change in moulin radius due to melting (m)
%uw  = water velocity at each node (m/s)
%Vadd = volume of water added due to melting for each grid cell (m3)


ks = 1; % meters, this is the mean hight of ice surface roughness elements...
% currently this is unconstrained...

include_ice_temperature = true; %true means that the change in the ice temperature is included in...
%the calculated change in moulin radius. If false, it makes the implicit
%assumption that the ice temperature and water temperature are both at the pressure melting temperature. 

Bathurst = true; %true means that the friction factor is calculated using..
%the Bathurst formulation... this equation is valid when
%roughness height ./ hydrualic diameter >= 0.05
% if false, the Colebrook-White formulation will be applied, which is only
% valid when roughness height ./ hydrualic diameter < 0.05


if include_ice_temperature
    Ti = importTz('HarrS4C', z);
else
    Ti = NaN; %#ok<UNRCH>
end


%% Initiate calculations for turbulent melting


waterpresent = hw - z; %logical matrix to determine in what nodes water is present 
waterpresent(waterpresent<0) =0;
S = (pi * Mr.^2);
uw = Qout ./  S; %calculate the water velocity in each node 
uw(waterpresent ==0) =0; % if there is no water in a given cell,

% ------------------------------
% Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity.
%It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross
% section.
if uw > 9.3
    disp('big velocity !!!')
    disp('assigning terminal velocity of 9.3ms-1')
end

uw = min(uw,9.3);
% ------------------------------

Mp   = 2 .* pi .* Mr; % wetted/melting perimeter
Dh   = (4*(pi .* Mr.^2)) ./ Mp; %hydrualic diameter
Rh   = (pi.* Mr.^2) ./ Mp; % hydraulic radius
% just give a flag if the ks/dh ratio is less than 0.05
% if (ks/Dh) < 0.05
%     disp('ks/dh < 0.05')
%     disp('Should be using Colebrook-White')
% end

% calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
Pw   = C.rhow .* C.g .* waterpresent;
Tmw  =  273.16 - 9.8e-8.*(Pw - 611.73);  


% select the appropriate parameterization of DW friction factor
if Bathurst
    fR = 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %Bathurst parameterization for DW friction factor
else
    %fR = (1./(-2.*log10((ks/Dh)./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
    fR = 20;
end


% calculate the lengthscale over which head loss is determined...
% this is simply the length of the model grid cells unless they become
% non-uniform
dz          = nanmean(diff(z));

% calculate head loss following Munson 2005
hL  =  ((uw.^2) .* fR .* dz) ./(2 .* Dh .* C.g);


if include_ice_temperature
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dz)) ...
                ./ (2 .* pi .* Mr .* C.rhoi .* (C.cw .* (Tmw - Ti) + C.Lf)); 
    %This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
    % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
    % include the surrounding ice temperature 
    
else
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dz)) ./ (2 .* pi .* Mr .* C.rhoi .* C.Lf); %#ok<UNRCH>
    % This parameterization is closer to that which is traditionally used
    % to calculate melting within a subglacial channel where the ice and
    % water are at the same temperature
end

% Make sure that there is no melting in places without water
dM_dt(waterpresent == 0) = 0;

dM  = dM_dt .* dt; %change in radius over the given time step

Vadd = C.rhoi/C.rhow * trapz(2*pi*Mr.*dM, z); %volume of meltwater gained due to melting the surrounding ice





end







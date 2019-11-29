function [dM, uw, Vadd, head_loss] =turbulence_hutter(hw, Qout, Mr, z, dt, C)

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


ks = 0.01; % meters, this is the mean hight of ice surface roughness elements...
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
    Ti = importTz('Luthi', z);
else
    Ti = NaN;
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
if (ks/Dh) < 0.05
    disp('ks/dh < 0.05')
    disp('Should be using Colebrook-White')
end

% calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
Pw   = C.rhow .* C.g .* waterpresent;
Tmw  =  273.16 - 0.098.*(Pw - 611.73);

% select the appropriate parameterization of DW friction factor
if Bathurst
    fR = 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %Bathurst parameterization for DW friction factor
else
    fR = (1/(-2.*log10((ks/Dh)./3.7))).^2; %Colebrook-White parameterization for DW friction factor
end



if include_ice_temperature
    dz          = nanmean(diff(z)); %find the length of the nodes to use as the length scale
    melt        = (Mp .* C.kw .* Nu .* (Tw - Ti) ./ (4 .* C.Lf .* Rh)) ./( dz);    
    dTwdz       = (diff(Tw) ./ diff(z));
    dTwdz       = [0, dTw];
    dTw         =  -uw .* dTwdz  + 1 ./ (C.rhow .* C.cp .* S) ...
            .* (Mp .* tau0 .* uw - mdot .* (C.Lf + C.cp .* ...
            (Tw - Ti) - ((uw.^2)./2)));
    
else
    dz          = nanmean(diff(z)); %find the length of the nodes to use as the length scale
    head_loss   = ((uw.^2) .* C.f_moulin .* dz) ./(2.* Dh .* C.g);
    melt        = (Qout .* C.rhow .* C.g .* (head_loss ./ dz))./ ...
             (2 .* pi .* C.Lf .* Mr); % these units are m per second of wall melt

    
end

dM  = melt .* dt; %change in radius over the given time step
Vadd = C.rhoi/C.rhow * trapz(2*pi*Mr.*dM, z); %volume of meltwater for each node

if max(dM) > 10
    disp('big melt error')
end
%
% Water volume change:
% Vadd = (C.rhoi ./ C.rhow) .*(pi .* ((Mr+dM) - Mr)); %volume of meltwater for each node
% Vadd = trapz(Vadd,z); % return a single number, not a vector
%
% KP update: the function needs to return a volume, Vadd, that is the
% volume of water added to the moulin's water.  It is a 3D volume.
%Vadd = C.rhoi/C.rhow * trapz(2*pi*Mr.*dM, z);


end







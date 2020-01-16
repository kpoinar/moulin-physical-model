function [dOC, uw_oc, Voc] =openchannel(hw, Qin, Mr, Mxu, dt, Ti, z, relative_roughness, Bathurst, include_ice_temperature,dz)


% so the hydrualic radius will be calculated as the perimeter of 1/2 an
% egg, where the minor radius is determined purely from dC+dE and the
% major radius is calculated as dC+dE+dOC.
% but that means the dM now has a differnt hydraulic radius for where the
%water comes and goes, so the top of the moulin is more oval, and the
%bottom of the moulin is more circular
%
% Moulin discharge (Qout) is dependent on the height of the water within the
% moulin (hw) and the capacity of the subglacial system to accept water.

%inputs from the moulin model
%hw   = height of the water in the moulin (m)
%Qout = discharge of water into the subglacial system (m3/s)
%Mr   = current moulin radius (m)
%Mxu  = position of the upstream boundary of the moulin wall
%dt   = timestep 
%Ti   = ice temperature 
%z    = elevation of each model node 0 = bottom of moulin; H= z(top of
        %moulin)
%relative_roughness = 
        %relative roughness needed to calculate head loss
%Bathurst =
        %choose the type of turbulence parameteriztion used


C = makeConstants;
ks = relative_roughness;

%exports
%dM  = change in moulin radius due to melting (m)
%uw  = water velocity at each node (m/s)
%Vadd = volume of water added due to melting for each grid cell (m3)



%% Initiate calculations for turbulent melting

% calculate the lengthscale over which head loss is determined...
% this is simply the length of the model grid cells unless they become
% non-uniform
 %This is the change in the vertical


% calculate the length overwhich the water is experiencing headloss
dX = diff(Mxu);
dX(end+1) =  dX(end);

%short length 
%dL = dX;

%long length;
dL = sqrt(dz.^2 + dX.^2);

waterpresent = hw - z; %logical matrix to determine in what nodes water is present 
waterpresent(waterpresent>=0) =0; %water is present
waterpresent(waterpresent<0) =1; %water is not present, so calculations need to be done.


S = (pi * Mr.^2) /2; % The cross-sectional area is 
uw = Qin ./  S; %calculate the water velocity in each node 
uw(waterpresent ==0) =0; % if there is water in the given cell we dont do this calculation

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

Mp   =  pi .* Mr; % wetted/melting perimeter
Dh   = (2*(pi .* Mr.^2)) ./ Mp; %hydrualic diameter
Rh   = (0.5*pi.* Mr.^2) ./ Mp; % hydraulic radius
    % flag if the ks/dh ratio is less than 0.05
    % if (ks/Dh) < 0.05
    %     disp('ks/dh < 0.05')
    %     disp('Should be using Colebrook-White')
    % end

% ------------------------------    
% calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
Pw   = C.rhow .* C.g .* waterpresent;
Tmw  =  273.16 - 9.8e-8.*(Pw - 611.73);  


% select the appropriate parameterization of DW friction factor
if Bathurst
    fR =10* 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %#ok<UNRCH> %Bathurst parameterization for DW friction factor
else
    fR = (1./(-2.*log10((ks./Dh)./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
    %fR(:,ii) = 0.01;
end




% calculate head loss following Munson 2005
hL  =  ((uw.^2) .* fR .* dL) ./(2 .* Dh .* C.g);


if include_ice_temperature
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dL)) ...
                ./ (2 .* pi .* Mr .* C.rhoi .* (C.cw .* (Tmw - Ti) + C.Lf)); 
    %This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
    % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
    % include the surrounding ice temperature 
    
else
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dL)) ./ (2 .* pi .* Mr .* C.rhoi .* C.Lf); %#ok<UNRCH>
    % This parameterization is closer to that which is traditionally used
    % to calculate melting within a subglacial channel where the ice and
    % water are at the same temperature
end

% Make sure that there is no melting in places without water
dM_dt(waterpresent == 0) = 0;

dM  = dM_dt .* dt; %change in radius over the given time step

Vadd = C.rhoi/C.rhow * trapz(2*pi*Mr.*dM, z); %volume of meltwater gained due to melting the surrounding ice





end







function [dM, uw, Vadd] =turbulence(hw, Qout, Mrminor, Mrmajor, Mxd, Ms, dt, Ti, dz, z, wet, relative_roughness, Bathurst, include_ice_temperature)


%
% Moulin discharge (Qout) is dependent on the height of the water within the
% moulin (hw) and the capacity of the subglacial system to accept water.

%inputs from the moulin model
%hw   = height of the water in the moulin (m)
%Qout = discharge of water into the subglacial system (m3/s)
%Ms   = current moulin cross-section area (m^2) replace Mr (changed by CT)
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


% waterpresent = hw - z; %logical matrix to determine in what nodes water is present 
% waterpresent(waterpresent<0) =0;
%S = (pi .* Mrminor .*Mrmajor); (CT) replace this by an S_m input for the
%moulin cross-section area
uw = Qout ./  Ms; %calculate the water velocity in each node 
% uw(waterpresent ==0) =0; % if there is no water in a given cell,
uw(~wet) = 0; % if there is no water in a given cell, there is no water velocity

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

Mp   = pi.* (3 .*(Mrminor + Mrmajor) - sqrt((3.* Mrminor + Mrmajor) .* (Mrminor +3 .* Mrmajor))); % wetted/melting perimeter =  ellipse perimeter approx pi [ 3(Mrminor+Mrmajor) - sqrt((3*Mrminor+Mrmajor)(Mrminor+3*Mrmajor))]
Dh   = (4.*(pi .* Mrminor .* Mrmajor)) ./ Mp; %hydrualic diameter
Rh   = (pi.* Mrminor .* Mrmajor) ./ Mp; % hydraulic radius
    % flag if the ks/dh ratio is less than 0.05
    % if (ks/Dh) < 0.05
    %     disp('ks/dh < 0.05')
    %     disp('Should be using Colebrook-White')
    % end

% ------------------------------    
% calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
Pw   = C.rhow .* C.g .* wet;
Tmw  =  273.16 - 9.8e-8.*(Pw - 611.73);  


% select the appropriate parameterization of DW friction factor
if Bathurst
    fR =10* 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %#ok<UNRCH> %Bathurst parameterization for DW friction factor
else
    fR = (1./(-2.*log10((ks./Dh)./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
    %fR(:,ii) = 0.01;
end

fR =1;

% calculate the lengthscale over which head loss is determined...
% this is simply the length of the model grid cells unless they become
% non-uniform
%%%%%%%%%%%%% calculate the length overwhich the water is experiencing headloss
 dL = diff(Mxd);
 dL = [dL(1); dL];
 dL(end) =  dL(end-1);
 dL = max(dL,0);  % protect against negative dL
 %dL = dz;% 
 dL = sqrt(dL.^2 + dz.^2);
%dz          = nanmean(diff(z));

% calculate head loss following Munson 2005
hL  =  ((uw.^2) .* fR .* dL) ./(2 .* Dh .* C.g);


if include_ice_temperature
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dL)) ...
                ./ (Mp .* C.rhoi .* (C.cw .* (Tmw - Ti) + C.Lf)); 
    %This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
    % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
    % include the surrounding ice temperature 
    
else
    dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dL)) ./ (Mp.* C.rhoi .* C.Lf); %#ok<UNRCH>
    % This parameterization is closer to that which is traditionally used
    % to calculate melting within a subglacial channel where the ice and
    % water are at the same temperature
end

% Make sure that there is no melting in places without water
dM_dt(~wet) = 0;

dM  = dM_dt .* dt; %change in radius over the given time step

% Smooth it over a 10 meter reach
nsm = round(10/ 1);
dM = fastsmooth(dM,nsm,3,1);

Vadd = C.rhoi/C.rhow * trapz(z, Mp .* dM) / dt; %volume of meltwater gained due to melting the surrounding ice





end







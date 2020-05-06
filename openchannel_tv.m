function [dOC, Vadd_oc, hL, dL, fR, neg2, Rh, cres_area, Mp] = openchannel_tv(hw, Qin, Mr_minor, Mr_major, Mxu, dt, Ti, dz, z, wet, include_ice_temperature,   fR_oc_variable, relative_roughness_OC,  fR_oc_fixed)


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


%exports
%dM  = change in moulin radius due to melting (m)
%uw  = water velocity at each node (m/s)
%Vadd = volume of water added due to melting for each grid cell (m3)



%% Initiate calculations for turbulent melting

% calculate the lengthscale over which head loss is determined...
% this is simply the length of the model grid cells unless they become
% non-uniform
 %This is the change in the vertical


%%%%%%%%%%%%%% zero the values of water below the water line, now done by
%%%%%%%%%%%%%% wet
% waterpresent = hw - z; %logical matrix to determine in what nodes water is present 
% waterpresent(waterpresent>=0) =0; %water is present
% waterpresent(waterpresent<0) =1; %water is not present, so calculations need to be done.

%%%%%%%%%%%%% calculate the length overwhich the water is experiencing headloss
dL = diff(Mxu);
 
 dL = [dL(1); dL];
 dL(end) =  dL(end-1);
 neg = dL;

 neg2 = neg;
neg(neg<0) = 0;
neg(neg>=0)=1;
 
 dL = max(dL,0);  % protect against negative dL - 
 
 %repl = dL<=0;
 
 %dL = dz;% 
 dL = sqrt(dL.^2 + dz.^2);

%%%%%% make the assumption that the terminal velocity is reached
%%%%%% immediately

term_vel = 9; 
 
% from the terminal velocity and Qin calculate the area of a cresent 
cres_area = Qin ./ term_vel; 
     
% wetted perimeter is still 1/2 that of an ellipse
[Mp, ~] = ellipseperimeter(Mr_minor, Mr_major);
Mp = Mp/2;




Dh   = (4 .* cres_area) ./ Mp; %hydraulic diameter
Rh   = cres_area ./ Mp; % hydraulic radius

%%%%%%%%%%%%% The water temperature is assumed to be zero because it is not,
% technically, under pressure of the water above it
Tmw  =  273.15;   

%%%%%%%%%%%%%% select the appropriate parameterization of DW friction
%%%%%%%%%%%%%% factor 
 if fR_oc_variable
    ks = relative_roughness_OC;
    %fR = (1./(-2.*log10((ks./Dh)./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
        %fR =100* 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %#ok<UNRCH> %Bathurst parameterization for DW friction factor
        %fR = 10*(1./ (-1.987 .* log(ks./(5.15.*Dh)))).^2; 
        fR = 1./(-2*log10(ks./(5.15*Rh))).^2;
 else 
     fR = fR_oc_fixed;
     
 end



%These should now be set in call_moulin_geom 
% % fR = 0.1;  % typical value of roughness used in subg models and suchlike
% %fR = 100;  % ~75 is maximum observation reported by Gulley et al. (2013)
% fR = 1;
% fR = 0.1;
% fR = 1e-9;



%%%%%%%%%%%%%% expected headloss based on the discharge

hL = (((term_vel).^2) .*fR .* dL) ./ (2 .* Dh .* C.g);

%%%%%%%%%%%%% calculate the expected melt 
if include_ice_temperature
    dOC_dt =    (C.rhow .* C.g .* Qin .* (hL./dL)) ...
                ./ (Mp .* C.rhoi .* (C.cw .* (Tmw - Ti) + C.Lf));

    %This is modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
    % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
    % include the surrounding ice temperature 
    
else
    dOC_dt =    (C.rhow .* C.g .* Qin .* (hL./dL)) ./ (Mp .* C.rhoi .* C.Lf); %#ok<UNRCH>
    % This parameterization is closer to that which is traditionally used
    % to calculate melting within a subglacial channel where the ice and
    % water are at the same temperature
end

% Make sure that there is no melting in places without water
dOC_dt(wet) = 0;

dOC  = dOC_dt .* dt;%change in radius over the given time step



%OC fix 1 (lca 4/23):  Smooth it over a 10 meter reach
% nsm = round(10/ 1);
% dOC = fastsmooth(dOC,nsm,3,1);
% 


% % Replace anywhere with dL=0 with the potential drop:
%   dP = 0.5 * C.rhow / C.rhoi *C.g / C.Lf * Qin*dt ./ (2*pi*Mr_major);
  
%this should zero out melting if there is a negative slope
dOC = dOC.*neg;




Vadd_oc = C.rhoi/C.rhow * trapz(z, Mp .* dOC); %volume of meltwater gained due to melting the surrounding ice

% hL_L = hL./dL;
% hL_L(waterpresent == 0) = 0;


end







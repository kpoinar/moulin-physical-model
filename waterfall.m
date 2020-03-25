function [dOC, Vadd_oc] = waterfall(Mxu       %hw, Qin, Mr_minor, Mr_major, Mxu, dt, Ti, dz, z, relative_roughness, Bathurst, include_ice_temperature, wet, )



%steps


% calculate the slope of the of the upstream sidewall
 dL = diff(M.xu);
 dL = [dL(1); dL];
 dL(end) =  dL(end-1);
 dL = sqrt(dL.^2 + dz.^2); %this is the slope of the 


%% define a set of if statements to call various functions


if dL 


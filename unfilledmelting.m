function [dOC, Vadd_oc] = waterfall(Mxu       %hw, Qin, Mr_minor, Mr_major, Mxu, dt, Ti, dz, z, relative_roughness, Bathurst, include_ice_temperature, wet, )



%steps


% calculate the slope of the of the upstream sidewall
 dL = diff(M.xu);
 dL = [dL(1); dL];
 dL(end) =  dL(end-1);
 dL = sqrt(dL.^2 + dz.^2); %this is the slope of the 


%% define a set of if statements to call various functions


if dL >= 0 
    
    % call openchannel
    
elseif dL < 0 
    % call potential drop
    
end


% now determine if there is a change from negative/zero to positive

test_mxu  = [1.8 1 1.2 1.5 1.8 1.9 1.3 1 0.8 1.8]';
test_dz   = [1 1 1 1 1 1 1 1 1 1]';
test_dL   = diff(test_mxu);
test_dL  = [ test_dL; test_dL(end)];
dL_sign  = sign(test_dL);
dL_sign(dL_sign==0) =1;
test_dL  = dL_sign .* sqrt(test_dL.^2 + test_dz.^2)



for ii = 1:length(test_dL)
    if test_dL(ii) >= 0
        dT(ii) = 2;
        
    elseif test_dL(ii) < 0
        dT(ii) = 1;
    end
end
   
test_1 = diff(dL_sign)

% need to do 
   
figure 
hold on
plot( test_mxu, [10:-1:1], 'r')

yyaxis right 
plot(test_dL, [1:10])



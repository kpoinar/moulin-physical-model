% function [dOC, Vadd_oc] = waterfall(Mxu       %hw, Qin, Mr_minor, Mr_major, Mxu, dt, Ti, dz, z, relative_roughness, Bathurst, include_ice_temperature, wet, )
% 
% 
% 
% %steps
% 
% 
% % calculate the slope of the of the upstream sidewall
%  dL = diff(M.xu);
%  dL = [dL(1); dL];
%  dL(end) =  dL(end-1);
%  dL = sqrt(dL.^2 + dz.^2); %this is the slope of the 
% 
% 
% %% define a set of if statements to call various functions
% 
% 
% if dL >= 0 
%     
%     % call openchannel
%     
% elseif dL < 0 
%     % call potential drop
%     
% end
% 
% %%
% % now determine if there is a change from negative/zero to positive
%             %bottom of moulin                %top of moulin
% clear tmpcheck
% close all
%             test_mxu  = [0 0.5 1.2 1.5 1.8 1.9 1.3 1 0.8 1.8]';
% test_dz   = [1 1 1 1 1 1 1 1 1 1]';
% height1 = 9:-1:0
% test_dL   = diff(test_mxu);
% test_dL  = [ test_dL; test_dL(end)];
% dL_sign  = sign(test_dL);
% dL_sign(dL_sign==0) =1;
% test_dL  = dL_sign .* sqrt(test_dL.^2 + test_dz.^2)
% 
% 
% tmp2 = nan(10,10);
% tmp_dz = nan(10,10);
% for ii = 1:10
%     tmp2(ii,1:end) = (  test_mxu(1:end) -test_mxu(ii))  ;
%     tmp_dz(ii,1:end) =  (  height1(1:end)- height1(ii) ) ;
% end
% 
% 
% for ii = 1:length(test_dL)
%     if test_dL(ii) >= 0
%         tmpcheck(ii) = 1; %apply potential drop
%         
%     elseif test_dL(ii) < 0   
%        for jj = 1:10   
%             if tmp2(ii,jj) >0 && tmp_dz(ii,jj) < -1
%                 tmpcheck(ii) = 2; %apply potential drop 
%             elseif tmp2(ii,jj) <0 && tmp_dz(ii,jj) > -1
%                 tmpcheck(ii) = 3; %apply waterfall
%             end
%        end
%     end
% end
%    
%    
% figure 
% hold on
% plot( test_mxu, [0:9], 'r')
% ylim([-1 10])
% yyaxis right 
% plot(tmpcheck, [0:9], 'o')
% ylim([-1 10])


%%

clear all
clear tmpcheck
close all
            test_mxu  = [0 0.5 1.2 2.1 2 1.8 1.3 1 0.8 1.8]';
test_dz   = [1 1 1 1 1 1 1 1 1 1]';
height1 = 9:-1:0;
test_dL   = diff(test_mxu);
test_dL  = [ test_dL; test_dL(end)];
dL_sign  = sign(test_dL);
dL_sign(dL_sign==0) =1;
test_dL  = dL_sign .* sqrt(test_dL.^2 + test_dz.^2);


tmp2 = nan(10,10);
tmp_dz = nan(10,10);
for ii = 1:10
    tmp2(ii,ii:end) = ( test_mxu(ii) -  test_mxu(ii:end) )  ;
    tmp_dz(ii,ii:end) =  ( height1(ii) - height1(ii:end) ) ;
end


%tmp2 = flipud(tmp2);
%tmp_dz = flipud(tmp_dz);
for ii = 1:length(test_dL)
    if test_dL(ii) >= 0
        tmpcheck(ii) = 1; %apply potential drop
        
    elseif test_dL(ii) < 0   
       for jj = 1:10   
            if tmp2(ii,jj) <0  
                tmpcheck(ii) = 2; %apply potential drop 
            elseif tmp2(ii,jj) >=0 && dL_sign(
                tmpcheck(ii) = 3; %apply waterfall
               
            elseif tmp2(ii,jj) >= 0
                tmpcheck(ii) = 4 %apply waterfall
                
            end
       end
    end
end
 

figure 
hold on
plot( test_mxu, [0:9], 'g')
plot(tmpcheck, [0:9], 'o')
ylim([-1 10])
yyaxis right 
plot(test_dL, [0:9], 'r')
ylim([-1 10])


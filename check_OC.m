
time = seriesresults{2};



r1 = time.M.r_minor;
r2 = time.M.r_major;



area_semicircle = (pi .* r1.^2)./2;
perim_semicircle = (2 * pi *r1)./2;

hrd_semicircle = area_semicircle ./perim_semicircle;

area_semiellipse = (pi .* r1 .* r2) ./2;
perim_semiellipse =  (pi.* (3 .*(r1 + r2) - sqrt((3.* r1 + r2) .* (r1 +3 .* r2))))./2;
%perim_semiellipse =  pi * sqrt((r1.^2 + r2.^2)/2);
hrd_semiellipse = area_semiellipse./perim_semiellipse;

%%
figure

subplot(1,5,1)

hold on 
plot(time.oc_fR(:,end), time.z)
title('fR')


subplot(1,5,2)
plot(time.oc_rh(:,end), time.z)
title('Rh')


subplot(1,5,3)
plot(time.oc_cres_area(:,end), time.z)
title('area used')

subplot(1,5,4)
plot(time.oc_Mp(:,end), time.z)
title('perimeter used')


%%
figure
hold on;
plot(time.M.r_major(:,end), time.oc_rh(:,end), '.')
xlabel('Mr major')
ylabel('Rh from ellipse')
%% figure out what is going on with hydraulic radius

figure
hold on
plot( hrd_semiellipse,time.z)
plot( hrd_semicircle,time.z)

%%
% figure; 
% hold on 
% plot(r1, time.z)
% plot(r2, time.z)
% yyaxis right
% 
% plot(time.oc_rh(:,end), time.z)

figure 
subplot(1,2,1)
hold on
plot(area_semiellipse(:,1:100:end), time.z)
plot(perim_semiellipse(:,1:100:end), time.z)


subplot(1,2,2)
hold on
plot(area_semicircle(:,1:100:end), time.z)
plot(perim_semicircle, time.z)



%%
figure
hold on
plot(area_semicircle(:,1000), perim_semicircle(:,1000))

plot(area_semiellipse(:,1:100:1000), perim_semiellipse(:,1:100:1000))

axis([0 4 0 35])
xlabel('area')
ylabel('perimeter')
legend('semicircle', 'semiellipse')

% figure 
% plot( r1)
% plot(area_semicircle, r2)

%%

% xx = 50
% figure
% subplot(1,3,1)
% hold on
% plot(time.dP(:,288*xx), time.z)
% plot([0 0.005],[time.hw(288*xx), time.hw(288*xx)] )
% yyaxis right 
% plot(time.wet(:,288*xx), time.z)
% axis([0 0.004, 0 700])
% %%
% 
% subplot(1,3,2)
% plot(time.dC_major(:,end), time.z)
% 
% 
% subplot(1,3,3)
% hold on
% %axis([0 0.015 0 700])
% plot([0 0.015], [488 488], 'r')
% plot([0 0.015], [time.hw(end) time.hw(end)], 'g')
% plot(time.dP(:,end), time.z, 'linewidth', 2)
% title('dP')




% create cosine functions for each elevation band 
clear variables

load ~/Documents/Data/Yang_basins/mean_runoff_elevationbands_new.mat

% elevation band for the prescribed moulins: are 7-800 (one); 8-900(19) ; 9-1000(52); 10-1100(70); 11-1200(92); 12-1300(117);
% 13-1400(112); 14-1500(153); 15-1600(122); 16-1700(20); 17-1800 (two)
mean_runoff = mean_runoff_elev;

mean_runoff(:,5) = mean_runoff(:,4);
mean_runoff(:,4) = mean_runoff(:,2) - mean_runoff(:,3);

figure
hold on 
for ii = 1:11
   plot([ii ,ii], [mean_runoff(ii,3), mean_runoff(ii,2)], 'linewidth',5 )
end

plot((1:11), mean_runoff(:,1), 'o', 'markerfacecolor', 'w', 'linewidth', 3, 'markeredgecolor', 'k')
axis([0 12 1 7])
text(0.93,mean_runoff(1,1), '1', 'fontsize', 20)
text(1.84,mean_runoff(2,1), '19', 'fontsize', 20)
text(2.84,mean_runoff(3,1), '52', 'fontsize', 20)
text(3.84,mean_runoff(4,1), '70', 'fontsize', 20)
text(4.84,mean_runoff(5,1), '92', 'fontsize', 20)
text(5.82,mean_runoff(6,1), '117', 'fontsize', 15)
text(6.82,mean_runoff(7,1), '112', 'fontsize', 15)

text(7.82,mean_runoff(8,1), '153', 'fontsize', 15)
text(8.82,mean_runoff(9,1), '122', 'fontsize', 15)
text(9.84,mean_runoff(10,1), '20', 'fontsize', 20)
text(10.93,mean_runoff(11,1), '2', 'fontsize', 20)
xlabel('Moulin elevation (km)', 'fontweight', 'bold')
ylabel('Qin (mean & spread; m^3 s^{-1})', 'fontweight', 'bold')
xticks(1:11)
xtickangle(45)
xticklabels({'0.7-0.8','0.8-0.9','0.9-1.0','1.0-1.1','1.1-1.2','1.2-1.3','1.3-1.4','1.4-1.5','1.5-1.6','1.6-1.7','1.7-1.8'})
%mean_runoff = round(mean_runoff.*10)


%% Now create the sine curves for each elevation band 

sec = 86400*365;
dt = 3600*24 *0.0625; % seconds (1 day here)
% tmax: time to run model
tmax = 2.3*sec;%1.7 * sec; % seconds (5 years here)
tmax = 1.5*sec;
time.t = dt:dt:tmax;



for ii = 1:11
    Q(ii,:) =(mean_runoff(ii,4)/2 )*cos((2*pi.*time.t./(86400)))+mean_runoff(ii,1);
end



figure
hold on 
plot(time.t/86400, Q, 'linewidth',2)
xlim([dt/86400 dt*100/86400])
xlabel('Time (days)', 'fontweight', 'bold')
ylabel('Qin  m^3 s^{-1})', 'fontweight', 'bold')


%%
x = 0:10000:500000
h = sqrt(((2*100000)/(9.81*910))*(x))

figure
plot(x/1000,h)
axis([0 100 0 1500])

%%

Q = Q';
Q(:,2:12) = Q;
Q(:,1) = time.t;
%%
% clear all
% cd ~/Documents/Repositories/moulin/moulin-physical-model/
% load foxxro_11_12.mat
% sec = 86400*365;
% dt = 3600*24 *0.0625; % seconds (1 day here)
% % tmax: time to run model
% tmax = 2.3*sec;%1.7 * sec; % seconds (5 years here)
% tmax = 1.5*sec;
% time.t = dt:dt:tmax;
% Qin = interp1((foxxro.date-foxxro.date(1))*86400,foxxro.runoff,time.t);
% Qin(isnan(Qin)) = 0;
% 
% 
% 
% 
% tmptime = 1:86400;
% 
% load foxxro_11_12.mat
% foxxtemp(:,1) = foxxro.date(3759:6190);
% foxxtemp(:,2) = foxxro.runoff(3759:6190);
% 
% Qfox = interp1((foxxtemp(:,1)-foxxtemp(1,1))*86400,foxxtemp(:,2),time.t);
% 
% %Qsin = 5* cosd(tmptime);
% Qsin = 2.*cos((2*pi.*time.t./(86400)))+3;
% Qin2 = Qsin + cos((2*pi.*time.t./(86400*5)))+1;
% Qin3 = 2.*cos((2*pi.*time.t./(86400))) + 3*cos((pi.*time.t./(86400*3)))+6;
% 
% figure
% hold on
% plot(time.t(1:48*4), Qsin(1:48*4))
% % plot(time.t(1:48*4), Qin2(1:48*4))
%  plot(time.t(1:48*4), Qin3(1:48*4))
% % plot(time.t(1:48*4), Qfox(1:48*4))
% 
% 
% 
% Qin4 =0.8*Qfox+3.5;
% %%
% 
% 
% clear tmp
% tmp = 1:-0.001:0.01;
% %tmp(101:249) = tmp(100);
% Qin4(1120:1623) = Qin4(1120:1623).*tmp(1:504);
% figure;
% plot(Qin4)
% 
% 
% Qcos2(:,1) = time.t;
% Qcos2(:,2) = Qsin;
% Qcos2(:,3) = Qin2;
% Qcos2(:,4) = Qin3;
% Qcos2(:,5) = Qfox;
% Qcos2(:,6) = Qin4;
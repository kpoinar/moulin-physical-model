clear all
cd ~/Documents/Repositories/moulin/moulin-physical-model/
load foxxro_11_12.mat
sec = 86400*365;
dt = 3600*24 *0.0625; % seconds (1 day here)
% tmax: time to run model
tmax = 2.3*sec;%1.7 * sec; % seconds (5 years here)
tmax = 1.5*sec;
time.t = dt:dt:tmax;
Qin = interp1((foxxro.date-foxxro.date(1))*86400,foxxro.runoff,time.t);
Qin(isnan(Qin)) = 0;




tmptime = 1:86400;

load foxxro_11_12.mat
foxxtemp(:,1) = foxxro.date(3759:6190);
foxxtemp(:,2) = foxxro.runoff(3759:6190);

Qfox = interp1((foxxtemp(:,1)-foxxtemp(1,1))*86400,foxxtemp(:,2),time.t);

%Qsin = 5* cosd(tmptime);
Qsin = 2.*cos((2*pi.*time.t./(86400)))+3;
Qin2 = Qsin + cos((2*pi.*time.t./(86400*5)))+1;
Qin3 = 2.*cos((2*pi.*time.t./(86400))) + 3*cos((pi.*time.t./(86400*3)))+6;

figure
hold on
plot(time.t(1:48*4), Qsin(1:48*4))
% plot(time.t(1:48*4), Qin2(1:48*4))
 plot(time.t(1:48*4), Qin3(1:48*4))
% plot(time.t(1:48*4), Qfox(1:48*4))



Qin4 =0.8*Qfox+3.5;
%%


clear tmp
tmp = 1:-0.001:0.01;
%tmp(101:249) = tmp(100);
Qin4(1120:1623) = Qin4(1120:1623).*tmp(1:504);
figure;
plot(Qin4)


Qcos2(:,1) = time.t;
Qcos2(:,2) = Qsin;
Qcos2(:,3) = Qin2;
Qcos2(:,4) = Qin3;
Qcos2(:,5) = Qfox;
Qcos2(:,6) = Qin4;
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



%Qsin = 5* cosd(tmptime);
Qsin = 2.*cos((2*pi.*time.t./(86400)))+3;
plot(time.t(1:48*4), Qsin(1:48*4))





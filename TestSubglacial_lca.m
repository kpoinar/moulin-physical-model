
clear variables 
%close all

%load MrAndZ_example.mat
H = 1000; %Ice thickness 
dz = 2;
z = 0:dz:H;
z = z';
Mr = 0*z +3; %ones(length(z),1)*5;

secinday = 3600*24;
dt = secinday ./24 ./ 4 ; % 5min
numdays = 50;
C = makeConstants;

%Mr_z = interp1(z,Mr,199);

opt = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));


H = 980; %Ice thickness
L = 15e3; %
y0 = [300 0.3];
%tspan = [0 secinday.*(1/24)]; % make the timespan1 hour
time1 = 0:dt:(secinday * numdays);
%Qin = ones(1,length(time1)) .* 3; %Recharge
load Qsine.mat
Qsine  = Qsine(8:end,:) ;
Qin         = interp1(Qsine(:,1), Qsine(:,2), time1, 'linear', 'extrap'); % run an interp just in case the timeframe changes
Qin(Qin<1) =1;
Qin =Qin*0.5 +1.5;


S = zeros(1,length(time1));
S(1) = 0.4;
hw(1) = 950;
Qout(1) = Qin(1);

for ii = 1:length(time1)
    
    tspan = [time1(ii) time1(ii)+dt];
    y0    = [hw(ii) S(ii)];
     looptest = ii;
    [hw(ii+1),S(ii+1),Qout(ii+1)] = subglacialsc(Mr,z,Qin(ii),H,L,C,tspan,y0,opt);

    if hw(ii+1) > H
        hw(ii+1) = H;
        
        
    elseif isnan(hw(ii+1))
        hw(ii+1) = H;
        
        Qout(ii+1)  = C.c3 .* (S(ii).^(5/4)) .* (((C.rhow .* C.g .* H) ./L).^ (0.5)) ... %Schoof SI Equation 6
                        ; %some estimate of sheet opening (?) 
        S(ii+1)  = S(ii) + C.c1 .* Qout(ii) .* ((C.rhow .* C.g .* H) ./L)...
                   - C.c2 .* ((C.rhoi .* C.g .* H - C.rhow .* C.g .* H).^C.n) .* S(ii); % Schoof SI equation 1 without sliding opening
    else
        hw(ii+1) = hw(ii+1);
        
    end
     
     if Qout(ii+1) < 0
         Qout(ii+1) = 0;
         ii
     end
end
% [t,y] = ode23s(@(t,y) subglacial_odefcn(t,y,Mr,z,Qin,H,dx,C),tspan, y0);   
% hw = y(:,1); % moulin head (m)
% S = y(:,2); % channel cross-section area
% Qout = C.c3 .* S(end).^(5/4) .* sqrt( C.rhow*C.g*hw(end)/dx); % discharge out the channel

% subplot(1,2,1)
% plot(hw)
% hold on
% %plot(t/secinday,S)
% subplot(1,2,2)
% plot(Mr, z)
% hold on
% plot(-Mr, z)
%%

figure
hold on
plot(S)
ylabel('S')
axis([0 1000 0 14])
yyaxis right
plot(hw)

axis([0 1000 0 500])
ylabel('hw')

figure
hold on
plot(Qin)
plot(Qout)
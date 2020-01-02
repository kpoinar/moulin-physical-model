
%load MrAndZ_example.mat
H = 1000; %Ice thickness 
dz = 2;
z = 0:dz:H;
z = z';
Mr = 0*z +3; %ones(length(z),1)*5;

secinday = 3600*24;
C = makeConstants;

%Mr_z = interp1(z,Mr,199);

Qin = 3; %Recharge
H = 1000; %Ice thickness
L = 15e3; %
y0 = [300 0.3];
tspan = [0 secinday.*(1/24)]; % make the timespan1 hour


[hw,S,Qout] = subglacialsc(Mr,z,Qin,H,L,C,tspan,y0)

load Qcosines.mat %1 = time, 2 cosine function, 3 
Qcos2   = Qcos2(1:end,:) ;
Qin = Qcos2(:,2);

%%
numofdays = 3;
ec       = 86400*numofdays;   %seconds * days
dt        = 900;        % Timestep, seconds (15 minutes)
tmax      = sec;        % seconds (5 years here)
time.t    = dt:dt:tmax; % seconds

load Qcosines.mat %1 = time, 2 cosine function, 3 
Qcos2   = Qcos2(1:end,:) ;
%change Qcos2 column for different types: 1. cosine, 2. cosine with small
%melt event, 3. cosine with large melt event, 4. quasi-real data
Qin     = interp1(Qcos2(:,1), Qcos2(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
Qin     = Qin*0.8 +3; %scale Qin to deal with a few model issues



for ii = 2:length(Qin)
    tspan = [time.t(ii)
    [hw(ii),S(ii),Qout(ii)] = subglacialsc(Mr,z,Qin(ii),H,L,C,tspan,y0);
    y0 = [hw(ii-1), S(ii-1)];
    
end

% [t,y] = ode23s(@(t,y) subglacial_odefcn(t,y,Mr,z,Qin,H,dx,C),tspan, y0);   
% hw = y(:,1); % moulin head (m)
% S = y(:,2); % channel cross-section area
% Qout = C.c3 .* S(end).^(5/4) .* sqrt( C.rhow*C.g*hw(end)/dx); % discharge out the channel

subplot(1,2,1)
plot(hw)
hold on
%plot(t/secinday,S)
subplot(1,2,2)
plot(Mr, z)
hold on
plot(-Mr, z)
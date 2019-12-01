
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


[hw,S,Qout] = subglacialsc(Mr,z,Qin,H,L,C,tspan,y0);

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
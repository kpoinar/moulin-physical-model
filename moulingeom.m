% Solve for moulin geometry, considering
% 1. elastic deformation (closure)
% 2. creep deformation (closure)
% 3. refreezing (closure)
% 4. turbulent melting (opening)
% 5. ice lid?
% 6. xz-shear deformation?
%
%
clear variables
close all

set(0,'DefaultFigureWindowStyle','docked')

sec = 86400*365;
C = makeConstants;
Tdatatype = 'HarrS4C';
chebx = 0;  % chebx=1 is not working yet
nt = 1000;    % plot every nt timesteps
artesian = 1;  % allow moulin to shed water?

Qoutwinter = 0.01;  % minimum outflux (e.g. wintertime outflux)
% ^^ Just an estimate (m2/s)
Qscale = 1.0;     % factor to scale FOXX runoff by
ndaylag = 1/24;      % How many days to lag the Qin, Qout by?

E = 1.5;       % enhancement factor for creep
% E = 1;       % no enhancement




%
% H: the ice thickness
H = 700; % meters
R0 = 2;  % radius of moulin initially
% dt: timestep
dt = 3600*24 *0.0208333333333333; % seconds - now 0.5h
% tmax: time to run model
tmax = 2.3*sec;%1.7 * sec; % seconds (5 years here)
tmax = 0.5*sec; % a year and a half
% vertical spacing
dz = 1; % meters
% Minimum moulin width
Mrmin = 1e-3;  % 1 mm
% Prescribe an annual date of hydrofracture?  # if yes. Really high # if no.
HFdoy = 99999999999999;%165; % Mid June
%
% hw: height of water
hw = 1 * H * C.rhoi/C.rhow;  % Celia need to change this with subglacial model


% z: the vertical coordinate system, positive upward
z = (0:dz:H)';
nz = numel(z);
% t: time
time.t = dt:dt:tmax; % seconds
%
% Assign hydraulic gradient parameters and water flow speed:
ubottom = 1;  % m/s; this is something like 1 mm/sec
utop = ubottom; % just make something up for now
u = linspace(ubottom,utop,nz)';
u0 = ubottom; z0 = 0;
L = 100e3; % length scale over which to take hydraulic gradient
% dhdx = H / (30e3);    % just made something up
%
% Mr: the moulin radius vs depth (cylindrical symmetry)
Mr = R0*ones(size(z));
%Mr = linspace(0.1,R0,nz)';
%Mr = linspace(R0,0.1,nz)';
% R0(j) = 10*rand(1);
% Mr = fastsmooth(R0(j)*rand(size(z)),10,3,1);
Mrinit = Mr;
% How much water is in the moulin?
i = z<=hw;
V0 = trapz(z(i),pi*Mr(i).^2);
V = V0; Vturb = 0; Vfrz = 0;
% Pw: water pressure at flotation
Pw = C.rhow*C.g*hw;  
% Moulin water volume
% Mwv = pi*trapz(z,Mr.^2); % assume it is entirely filled
% Mwv = Mwv * hw/H; % now drop it to flotation
%
%
% Tz: temperature profile in ice
Tz = importTz(Tdatatype,z);
% Tair: air temperature timeseries
Tair = C.T0 - 8 - cos(2*pi*time.t/sec)*12; % Kelvin at every timestep
%Tair = C.To - 8 - cos(2*pi*(time.t+0.381*sec)/sec)*12;  % Start time at the melt season onset
% T: ambient ice temperature
Tfar = Tz; % Kelvin
%
% Set up x grid for far-field temperature
xmax = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx] = setupx(dt,chebx,xmax);
% 
T = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1) = C.T0;   % Melting point at the moulin wall
%
% Assign elastic deformation parameters
sigx = -50e3;%100e3;
sigy = -50e3;%-100e3;
tauxy = 100e3;%100e3;
%
%
%
figure(3); clf;
NP = 4; ii=1;
subplot(1,NP,1); hold on; title('Moulin radius'); xlabel('meters'); ii=ii+1;
subplot(1,NP,ii); hold on; title('Creep'); xlabel('m / dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Refreezing'); xlabel('m / dt'); ii=ii+1;
subplot(1,NP,ii); hold on; title('Turbulent melt'); xlabel('m /dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Potential energy melt-out'); xlabel('m /dt'); ii=ii+1;
subplot(1,NP,ii); hold on; title('Elastic'); xlabel('m / dt')
%
% Record minimum and maximum moulin diameter
time.Mrmm = zeros(2,numel(time.t));
time.V    = zeros(size(time.t));
time.hw   = zeros(size(time.t));
time.Vcapacity = zeros(size(time.t));
time.Vfrz = zeros(size(time.t));
time.Vturb = zeros(size(time.t));
dF = zeros(size(z));
dE = zeros(size(z));
dC = zeros(size(z));
dM = zeros(size(z));
dP = zeros(size(z));
%
%Qout = ubottom*pi*R0^2;
%const = (time.t(end)*Qout) / trapz(time.t,(Tair>C.To) .* (Tair-C.To));
% Expect Qin to be ~ 10 m3/s
const = 10 / (max(Tair-C.T0));
% Construct an approximate Qin for each timestep, based on air temps:
% Qin = double(Tair>C.To) .* (Tair-C.To)*const;
% Qin = zeros(size(time.t));


% real supraglacial discharge
%load foxxro_11_12.mat
%Qin = interp1((foxxro.date-foxxro.date(1))*86400,foxxro.runoff,time.t);
%Qin(isnan(Qin)) = 0;

%sinusoidal input
load Qsine.mat
Qin = interp1(Qsine(:,1), Qsine(:,2), time.t); % run an interp just in case the timeframe changes

% Scale the runoff up or down (usually down)
Qin = Qin * Qscale;



%
% Readjust Qout to be Qin but smoothed and lagged by some days
nlag = round(ndaylag*24*3600/dt);
Qout = fastsmooth(Qin,nlag/2,3,1);
% Qout = max(0,Qout);
Qout = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
Qout = Qout * 1;
% In the winter, there is still Qout.  Make it so.
Qout = max(Qout, Qoutwinter);
%Qout = zeros(size(time.t));
%
%
%
%
% Step through time
cc = 0;

for t = time.t
    cc = cc+1;
    % Use the PREVIOUS moulin radius in all calculations in each timestep,
    % so that the final result is not dependent on the order in which we do
    % creep, refreeze, turbulent melt, elastic, etc.
    Mrprev = Mr;
    %
    % Check if today is the day that we hydrofracture and reopen the bottom
    % of the moulin:
    if ~mod(HFdoy*86400 - t,sec)
        fprintf('Hydrofracture event! at t=%1.2f years (cc=%d)\n',t/sec,cc)
        Mr = max(Mr,R0);
    end
    %
    % Find the water level in the moulin at this timestep
    % Moulin water volume:
    V = watervolume(V,Vturb,Vfrz,Qin(cc),Qout(cc),dt);
    % Make an artesian spring if allowed and if V > Vmoulin
    if artesian, V = min(trapz(z,pi*Mr.^2),V); end
    % How high does that fill the moulin?
    hw = waterlevel(Mr,z,V);
            time.hw(cc) = hw;
            time.V(cc) = V;%trapz(V,z);
    hw_out(cc) = hw;        
    %
    %
    % Creep deformation: do this first because it is a larger term  
    % Change the water level twice-daily
    %     if mod(cc,0) % afternoon
    %         hw = 0.95*H; % above flotation
    %     else % at flotation
    %         hw = H * C.rhoi/C.rhow;  % meters
    %     end
    dC = creep(Mrprev,z,H,hw,T,dt,E,C);
    %
    % Refreezing
    %  T(z>hw,1) = Tair(cc);
    %  [~,dF,T,Vfrz] = refreeze(Mrprev,T,z,hw,dF,nx,x,dx,dt);
    %         time.Vfrz(cc) = Vfrz;
    %
    %
    % Turbulent melting: 
    % Water velocity at the bottom
    u0 = -Qout(cc) / (pi*Mr(1)^2*dz);
            time.u0(cc) = u0;
    % water velocity in the column
    %u = conserveWaterMass(Mr,z,u0,z0);
    
    % Turbulent melting
    Ti = C.T0 * ones(size(z)); % ice temperature: 0∞C
    Tw = C.T0 * ones(size(z)); % water temperature: 0∞C
    % Hutter Turbulence
    [dM(:,cc), u(:,cc), Vturb, head_loss(:,cc)] = turbulence_headloss(hw, Qout(cc), Mrprev, z, dt, C);  % dM(:,cc) for bug fixing, in standard run should be dM lca 11/18
            time.Vturb(cc) = Vturb;%trapz(Vturb,z);

    % Elastic deformation: do this last because it is a function of moulin 
    % radius.  Elastic deformation is small and sensitive to water pressure
    dE = elastic(z,Mrprev,hw,H,sigx,sigy,tauxy,C);

    %
    % Potential energy-based melting above the water line
    %dP = potentialdrop(Qin(cc),z,hw,Mrprev,dt,C);

    %
    % Now actually sum all the contributions to moulin size:
  %  Mr = Mr + dC + dF + dM(:,cc) + dE + dP; % dM(:,cc) for bug fixing, in standard run should be dM lca 11/18
  %  Mr = max(Mr,Mrmin);
  %  Mr_out(:,cc) = max(Mr,Mrmin); %lca bug fixing   
    Mr = Mr + dC + dF + dE + dP + dM(:,cc); %for bug fixing, in standard run should be dM lca 11/18
    
    Mr = max(Mr,Mrmin);
    Mr_out(:,cc) = max(Mr,Mrmin); %lca bug fixing 
    
    % Record moulin max and min radius at every timestep
    time.Mrmm(:,cc) = [min(Mr) max(Mr)];
    %
    % Record volume capacity of moulin
    time.Vcapacity(cc) = trapz(z,pi*Mr.^2);
    %
    % Plot profiles sometimes
    if ~mod(cc,nt)
        figure(3); ii=1;
            [~,jw] = min(abs(z-hw));
            subplot(1,NP,1); plot(Mr,z,Mr(jw),hw,'*k'); title(sprintf('Moulin radius at t=%1.1f days',t/3600/24)); ii=ii+1;
            subplot(1,NP,ii); plot(dC,z,dC(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dF,z,dF(jw),hw,'*k'); ii=ii+1;
            subplot(1,NP,ii); plot(dM,z,dM(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dP,z,dP(jw),hw,'*k'); ii=ii+1;
            subplot(1,NP,ii); plot(dE,z,dE(jw),hw,'*k');
    end
    %
end
%%
figure(30); clf; 
h1 = subplot(4,5,1:4);
plot(time.t/sec,time.Mrmm); legend('min Mr','max Mr')
ylabel('(m)')
title('Moulin diameter')

h2 = subplot(4,5,6:9);
plot(time.t/sec,time.hw);
title('Water level in moulin')
ylabel('(m)')
%set(gca,'yscale','log')

h3 = subplot(4,5,11:14);
plot(time.t/sec,Qin*dt,time.t/sec,Qout*dt);%,time.t/sec,time.V)
title(sprintf('Volume entering and leaving per dt (Q$_{in,out}$*dt), %1.2f day lag',ndaylag))
ylabel('(m$^3$)')
legend('Q$_{in}$ *dt','Q$_{out}$ *dt')%,'V in moulin')
%xlabel('Time (yrs)')

h4 = subplot(4,5,16:19);
plot(time.t/sec,time.Vturb,time.t/sec,-time.Vfrz); hold on
plot(time.t/sec,time.V,time.t/sec,time.Vcapacity)%'color',[0.75 0.75 0.3]);
ylabel('(m$^3$)')
xlabel('Time (yrs)')
legend('Vturb','-Vfrz','Vwater','Vcapacity')
title('Volumes')

linkaxes([h1, h2, h3, h4],'x');

subplot(4,5,[5 10 15 20])
[~,j] = min(abs(z-hw));
% bed
patch(2*max(Mr)*[-1 -1 1 1 -1],[-100 0 0 -100 -100],[0.9 0.7 0.4]); hold on
% water
patch([-Mr(1:j); flipud(Mr(1:j)); -Mr(1)],[z(1:j); flipud(z(1:j)); 0],[0.4 0.8 1]); hold on
set(gca,'xlim',2*max(Mr)*[-1 1],'ylim',[-100 H+100]);
% moulin walls
plot(Mr,z,'-k',-Mr,z,'-k')%,[-1 1]*Mr(j),[1 1]*hw,'-c')
% ice sheet surface
plot([-2*max(Mr) -Mr(end) NaN Mr(end) 2*max(Mr)],[1 1 1 1 1]*H,'-k')
xlabel('Moulin radius (m)')
ylabel('z (m)')
set(gca,'yaxislocation','right')
title('Final geometry')

%%


colors = brewermap(length(Qin), 'spectral');

figure
subplot(1,2,1)
addToolbarExplorationButtons(gcf)
hold on; box on; grid on 
for ii = 1:10:length(Qin)
plot( Mr_out(:,ii), (z), 'color', colors(ii,:))
end

subplot(1,2,2)
hold on; box on; grid on 
for ii = 1:10:length(Qin)
plot( u(:,ii), (z), 'color', colors(ii,:))
end

figure
hold on
plot(time.t, hw_out)

axis([time.t(1), time.t(24*5), 0 500])
%%
figure; 
plot(Qin); hold on; 
ylabel('Qin')
yyaxis right; 
%plot(Mr_out(250,:))
ylabel('moulin radius @ 250m')
axis([2000 4000 0 1400])

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
C = makeConstants;

% Model duration and time characterisitcs
sec = 86400*90;
tmax = sec .* 1 ; % 90 days %duration of the model run
dt = 3600*24 *(0.0208333333333333/2); % seconds - now 0.25h %TIMESTEP
time.t = dt:dt:tmax; % seconds


% Various tweakable parameters
Tdatatype = 'HarrS4C';
chebx     = 0;  % chebx=1 is not working yet
nt        = 1000;    % plot every nt timesteps
artesian  = 1;  % allow moulin to shed water?
HFdoy = 99999999999999;%165; % Mid June

ndaylag   = 1/24;      % How many days to lag the Qin, Qout by?


% Qin and Qout
load Qsine.mat
Qscale  = 1;
Qin     = interp1(Qsine(:,1), Qsine(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
Qin     = Qin * Qscale;

% Readjust Qout to be Qin but smoothed and lagged by some days
nlag    = round(ndaylag*24*3600/dt);
Qout    = fastsmooth(Qin,nlag/2,3,1);
Qout    = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
Qoutwinter ...
        = 0.01;  % minimum outflux (e.g. wintertime outflux) m3/s
Qout    = max(Qout, Qoutwinter);


%moulin characteristics
H       = 700; % meters % H: the ice thickness
R0      = 2;  % radius of moulin initially
Mrmin   = 1e-3;  % minimum moulin radius, meters

dz      = 1; % vertical spacing % meters
z       = (0:dz:H)'; %index 1 is the bottom of the moulin, index(H) is ice surface
E       = 1.0;       % enhancement factor for creep, 1.0 = no enhancement

%initializated moulin info
hw      = 1 * H * C.rhoi/C.rhow;  % Celia need to change this with subglacial model
Pw      = C.rhow*C.g*hw;  %Pw: water pressure at flotation
Mr      = ones(size(z)); %initialize moulin geometry 
Mrinit  = R0*Mr(:,1); %intial moulin geometry


% ice characteristics
Tz      = importTz(Tdatatype,z);
Tfar    = Tz; % Kelvin
xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx]...
        = setupx(dt,chebx,xmax);
T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1)  = C.T0;   % Melting point at the moulin wall


% Moulin water volume
i       = z<=hw;
V0      = trapz(z(i),pi*Mr(i).^2);
V       = V0; 
Vturb   = 0; 
Vfrz    = 0;
% 
Pw = C.rhow*C.g*hw;  

%  Elastic deformation parameters
sigx = -50e3;%100e3;
sigy = -50e3;%-100e3;
tauxy = 100e3;%100e3;
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

%% Run the model
% Step through time
cc = 0;

for t = time.t
    cc = cc+1;
    % Use the PREVIOUS moulin radius in all calculations in each timestep,
    % so that the final result is not dependent on the order in which we do
    % creep, refreeze, turbulent melt, elastic, etc.
    Mrprev = Mr;
 
%     % Check if today is the day that we hydrofracture and reopen the bottom
%     % of the moulin:
%     if ~mod(HFdoy*86400 - t,sec)
%         fprintf('Hydrofracture event! at t=%1.2f years (cc=%d)\n',t/sec,cc)
%         Mr = max(Mr,R0);
%     end
    
% Find the water level in the moulin at this timestep
    % Moulin water volume:
    V = watervolume(V,Vturb,Vfrz,Qin(cc),Qout(cc),dt);
    % Make an artesian spring if allowed and if V > Vmoulin
    if artesian, V = min(trapz(z,pi*Mr.^2),V); end
    % How high does that fill the moulin?
    [hw, cumvol, M0, min_index] = waterlevel(Mr,z,V);
       time.hw(cc)  = hw;
       time.V(cc)   = V;
      


% Creep deformation: do this first because it is a larger term  
    % Change the water level twice-daily
    dC = creep(Mrprev,z,H,hw,T,dt,E,C);
    time.dC(:,cc)     = dC; %preserve the creep change
    %
    % Refreezing
    %  T(z>hw,1) = Tair(cc);
    %  [~,dF,T,Vfrz] = refreeze(Mrprev,T,z,hw,dF,nx,x,dx,dt);
    %         time.Vfrz(cc) = Vfrz;
    %
    %
% Turbulent melting: 
    %u = conserveWaterMass(Mr,z,u0,z0);
    % physically based turbulence
    [dM, u, Vturb, head_loss] = turbulence_headloss(hw, Qout(cc), Mrprev, z, dt);  % dM(:,cc) for bug fixing, in standard run should be dM lca 11/18
            time.Vturb(cc) = Vturb;
            time.dM(:,cc)  = dM;
            time.uw(:,cc)  = u;
            time.head_loss(:,cc) ...
                           = head_loss;
    % Elastic deformation: do this last because it is a function of moulin 
    % radius.  Elastic deformation is small and sensitive to water pressure
    dE = elastic(z,Mrprev,hw,H,sigx,sigy,tauxy,C);
            time.dE(:,cc)   = dE;

    %
    % Potential energy-based melting above the water line
    %dP = potentialdrop(Qin(cc),z,hw,Mrprev,dt,C);

    %
    % Now actually sum all the contributions to moulin size:
     
    Mr = Mr + dC + dF + dE + dP + dM; %for bug fixing, in standard run should be dM lca 11/18
            time.Mr(:,cc) = Mr;
            
    Mr = max(Mr,Mrmin);
            time.Mr_minapplied(:,cc) = Mr;

    % Record moulin max and min radius at every timestep
            time.Mrmm(:,cc) = [min(Mr) max(Mr)];
    %
    % Record volume capacity of moulin
    time.Vcapacity(cc) = trapz(z,pi*Mr.^2);
    %

    %
end
%%
NPanels = 4; pan = 1;
spacing = 96; % daily with 15 min time steps
color1  = brewermap(length(time.t), '*spectral');


figure
subplot(1, NPanels,pan)
hold on
title('elastic')
xlabel('meters')
for ii = 1:spacing: length(time.t)
plot(time.dE(:,ii),z, 'color', color1(ii,:))

end




subplot(1, NPanels,pan+1)
hold on
title('Creep')
xlabel('meters')
for ii = 1:spacing: length(time.t)
plot(time.dC(:,ii),z, 'color', color1(ii,:))

end

subplot(1, NPanels,pan+2)
hold on
title('Turbulent melting')
xlabel('meters')
for ii = 1:spacing: length(time.t)
plot(time.dM(:,ii),z, 'color', color1(ii,:))

end

subplot(1, NPanels,pan+3)

hold on 
title('Moulin radius')
xlabel('meters')
ylabel('purple early, red late')
for ii = 1:spacing: length(time.t)
plot(time.Mr(:,ii),z, 'color', color1(ii,:))
end
%%
% 
% NP = 4; ii=1;
% subplot(1,NP,1); hold on; title('Moulin radius'); xlabel('meters'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Creep'); xlabel('m / dt'); ii=ii+1;
% % subplot(1,NP,ii); hold on; title('Refreezing'); xlabel('m / dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Turbulent melt'); xlabel('m /dt'); ii=ii+1;
% % subplot(1,NP,ii); hold on; title('Potential energy melt-out'); xlabel('m /dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Elastic'); xlabel('m / dt')
% 
% %%
% 
% 
%     % Plot profiles sometimes
%     if ~mod(cc,nt)
%         figure(3); ii=1;
%             [~,jw] = min(abs(z-hw));
%             subplot(1,NP,1); plot(Mr,z,Mr(jw),hw,'*k'); title(sprintf('Moulin radius at t=%1.1f days',t/3600/24)); ii=ii+1;
%             subplot(1,NP,ii); plot(dC,z,dC(jw),hw,'*k'); ii=ii+1;
% %             subplot(1,NP,ii); plot(dF,z,dF(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dM,z,dM(jw),hw,'*k'); ii=ii+1;
% %             subplot(1,NP,ii); plot(dP,z,dP(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dE,z,dE(jw),hw,'*k');
%     end
%     
% 
% %%
% 
% 
% 
% 
% figure(30); clf; 
% 
% h1 = subplot(4,5,1:4);
% plot(time.t/sec,time.Mrmm); legend('min Mr','max Mr')
% ylabel('(m)')
% title('Moulin diameter')
% 
% h2 = subplot(4,5,6:9);
% plot(time.t/sec,time.hw);
% title('Water level in moulin')
% ylabel('(m)')
% %set(gca,'yscale','log')
% 
% h3 = subplot(4,5,11:14);
% plot(time.t/sec,Qin*dt,time.t/sec,Qout*dt);%,time.t/sec,time.V)
% title(sprintf('Volume entering and leaving per dt (Q_{in,out}*dt), %1.2f day lag',ndaylag))
% ylabel('(m$^3$)')
% legend('Q_{in}*dt','Q_{out} *dt')%,'V in moulin')
% %xlabel('Time (yrs)')
% 
% h4 = subplot(4,5,16:19);
% plot(time.t/sec,time.Vturb,time.t/sec,-time.Vfrz); hold on
% plot(time.t/sec,time.V,time.t/sec,time.Vcapacity)%'color',[0.75 0.75 0.3]);
% ylabel('(m$^3$)')
% xlabel('Time (yrs)')
% legend('Vturb','-Vfrz','Vwater','Vcapacity')
% title('Volumes')
% 
% linkaxes([h1, h2, h3, h4],'x');
% 
% subplot(4,5,[5 10 15 20])
% [~,j] = min(abs(z-hw));
% % bed
% patch(2*max(Mr)*[-1 -1 1 1 -1],[-100 0 0 -100 -100],[0.9 0.7 0.4]); hold on
% % water
% patch([-Mr(1:j); flipud(Mr(1:j)); -Mr(1)],[z(1:j); flipud(z(1:j)); 0],[0.4 0.8 1]); hold on
% set(gca,'xlim',2*max(Mr)*[-1 1],'ylim',[-100 H+100]);
% % moulin walls
% plot(Mr,z,'-k',-Mr,z,'-k')%,[-1 1]*Mr(j),[1 1]*hw,'-c')
% % ice sheet surface
% plot([-2*max(Mr) -Mr(end) NaN Mr(end) 2*max(Mr)],[1 1 1 1 1]*H,'-k')
% xlabel('Moulin radius (m)')
% ylabel('z (m)')
% set(gca,'yaxislocation','right')
% title('Final geometry')
% 
% %%
% 
% 
% colors = brewermap(length(Qin), 'spectral');
% 
% figure
% subplot(1,2,1)
% addToolbarExplorationButtons(gcf)
% hold on; box on; grid on 
% for ii = 1:10:length(Qin)
% plot( Mr_out(:,ii), (z), 'color', colors(ii,:))
% end
% 
% subplot(1,2,2)
% hold on; box on; grid on 
% for ii = 1:10:length(Qin)
% plot( u(:,ii), (z), 'color', colors(ii,:))
% end
% 
% figure
% hold on
% plot(time.t, hw_out)
% 
% axis([time.t(1), time.t(24*5), 0 500])
% %%


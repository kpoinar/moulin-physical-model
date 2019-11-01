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
sec = 86400*365;
C = makeConstants;
Tdatatype = 'HarrS4C';
chebx = 0;  % chebx=1 is not working yet
nt = 1000;    % plot every nt timesteps
artesian = 0;  % allow moulin to shed water?





E = 5;       % enhancement factor for creep
E = 0.1;       % no enhancement




%
% H: the ice thickness
H = 1000; % meters
R0 = 2;  % radius of moulin initially
% dt: timestep
dt = 3600*24 *0.0625; % seconds (1 day here)
% tmax: time to run model
tmax = 2.3*sec;%1.7 * sec; % seconds (5 years here)
% vertical spacing
dz = 1; % meters
% Minimum borehole width
Mrmin = 1e-9;  % 1 mm
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
Tair = C.To - 8 - cos(2*pi*time.t/sec)*12; % Kelvin at every timestep
%Tair = C.To - 8 - cos(2*pi*(time.t+0.381*sec)/sec)*12;  % Start time at the melt season onset
% T: ambient ice temperature
Tfar = Tz; % Kelvin
%
% Set up x grid for far-field temperature
xmax = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx] = setupx(dt,chebx,xmax);
% 
T = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1) = C.To;   % Melting point at the moulin wall
%
% Assign elastic deformation parameters
sigx = -100e3;%100e3;
sigy = 0;%-100e3;
tauxy = 500e3;%100e3;
%
%
%
figure(3); clf;
subplot(1,6,1); hold on; title('Moulin radius'); xlabel('meters')
subplot(1,6,2); hold on; title('Creep'); xlabel('m / dt')
subplot(1,6,3); hold on; title('Refreezing'); xlabel('m / dt')
subplot(1,6,4); hold on; title('Turbulent melt'); xlabel('m /dt')
subplot(1,6,5); hold on; title('Potential energy melt-out'); xlabel('m /dt')
subplot(1,6,6); hold on; title('Elastic'); xlabel('m / dt')
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
const = 10 / (max(Tair-C.To));
% Construct an approximate Qin for each timestep, based on air temps:
% Qin = double(Tair>C.To) .* (Tair-C.To)*const;
% Qin = zeros(size(time.t));
load foxxro_11_12.mat
Qin = interp1((foxxro.date-foxxro.date(1))*86400,foxxro.runoff,time.t);
Qin(isnan(Qin)) = 0;



%Qin = Qin/100;



%
% Readjust Qout to be Qin but smoothed and lagged by some days
ndaylag = 18/24;
nlag = round(ndaylag*24*3600/dt);
Qout = fastsmooth(Qin,nlag/2,3,1);
% Qout = max(0,Qout);
Qout = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
Qout = Qout * 1;
%Qout = zeros(size(time.t));
%
%
%
%
% Step through time
cc = 0;

for t = time.t
    cc = cc+1;
    % Consider using the previous moulin radius in all calculations in each
    % timestep, so that the final result is not dependent on the order in
    % which I do creep, refreeze, turbulent melt, elastic, etc.
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
    % Make an artesian spring if V > Vmoulin
    if artesian, V = min(trapz(z,pi*Mr.^2),V); end
    % How high does that fill the moulin?
    hw = waterlevel(Mr,z,V);
            time.hw(cc) = hw;
            time.V(cc) = V;
    %
    
    % Creep deformation: do this first because it is a larger term  
    % Change the water level twice-daily
    %     if mod(cc,0) % afternoon
    %         hw = 0.95*H; % above flotation
    %     else % at flotation
    %         hw = H * C.rhoi/C.rhow;  % meters
    %     end
    dC = creep(Mrprev,z,H,hw,T,dt,E,C);
    
    % Refreezing
%     T(z>hw,1) = Tair(cc);
%     [~,dF,T,Vfrz] = refreeze(Mrprev,T,z,hw,dF,nx,x,dx,dt);
%             time.Vfrz(cc) = Vfrz;
    %
    
    % Turbulent melting: 
    % Water velocity at the bottom
    u0 = -Qout(cc) / (pi*Mr(1)^2*dz);
    time.u0(cc) = u0;
    % wate velocity in the column
    u = conserveWaterMass(Mr,z,u0,z0);
    % Turbulent melting
    [dM,Vturb] = turbulence(Mrprev,u,L,dt,H,hw,z,C);
            time.Vturb(cc) = Vturb;

    % Elastic deformation: do this last because it is a function of moulin 
    % radius.  Elastic deformation is small and sensitive to water pressure
    %dE = elastic(sigx,sigy,C.nu,C.E,tauxy,Pw,Mrprev);

    %
    % Potential energy-based melting above the water line
    %dP = potentialdrop(Qin(cc),z,hw,Mrprev,dt,C);

    %
    % Now actually sum all the contributions to moulin size:
    Mr = Mr + dC + dF + dM + dE + dP;
    Mr = max(Mr,Mrmin);
        
    % Record moulin max and min radius at every timestep
    time.Mrmm(:,cc) = [min(Mr) max(Mr)];
    %
    % Record volume capacity of moulin
    time.Vcapacity(cc) = trapz(z,pi*Mr.^2);
    %
    % Plot profiles sometimes
    if ~mod(cc,nt)
        figure(3); 
            [~,jw] = min(abs(z-hw));
            subplot(1,6,1); plot(Mr,z,Mr(jw),hw,'*k'); title(sprintf('Moulin radius at t=%1.1f days',t/3600/24))
            subplot(1,6,2); plot(dC,z,dC(jw),hw,'*k');
            subplot(1,6,3); plot(dF,z,dF(jw),hw,'*k'); 
            subplot(1,6,4); plot(dM,z,dM(jw),hw,'*k');
            subplot(1,6,5); plot(dP,z,dP(jw),hw,'*k');
            subplot(1,6,6); plot(dE,z,dE(jw),hw,'*k');
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
title(sprintf('Volume entering and leaving per dt (Q$_{in,out}$*dt), %1.1f day lag',ndaylag))
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
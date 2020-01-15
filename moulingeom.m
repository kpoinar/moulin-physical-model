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

warning('off', 'MATLAB:illConditionedMatrix')%These warnings off is not ideal, 
%but I havent figured out how to fix a problem with ode15s (couldnt find a solution for ode45, which just returned NaNs)
% - essentially when the water level is at or above the ice thickness,
% the initial conditions matrix becomes singlar/badly scaled - I *think*
% that this is okay in our circumstance (i.e. ode15s wants the hw to be
% greater than H, but it only converges at H, thought I am not totally sure
% on this) 
% The actual warning:
%Warning: Matrix is singular, close to singular or badly scaled. Results may be inaccurate. RCOND = NaN. 
%> In ode15s (line 589)
%  In subglacialsc (line 47)
%  In moulingeom (line 163)



% Do you want to plot using lauren's plotting function?
plot_using_lauren =  true;
save_figures      =  false;
visible_figures   =  true;

% Do you want to save the 'time' structure? This has (almost) everything
% documented...
save_timevariable =  false; 

save_location = '/Users/lcandre2/Documents/Projects/moulin_formation/modeloutputs';
datetime = datestr(now,'mm-dd-yyyy_HHMM'); %This will assign a unique date and time for both the figures and the model outputs
%% define som basic parameters
C         = makeConstants;  %constants used for parameterizations 
Tdatatype = 'Ryser_foxx';   %ice temperature profile to extrapolate from
numofdays = 10;             %set the number of days for the model run
H         = 500;            % ice thickness, meters
R0        = 3;              % radius of moulin initially
L         = 12e3;           % Length of the subglacial channel
f         = 0.5;            % fraction of the potential energy used to open the top of the moulin (above water level)

%inital guesses for subglacial model
hw(1) = H;                  % moulin water level (m)
S(1)  = R0;                 % subglacial channel cross sectional area (m^2)


chebx     = 0;              % chebx=1 is not working yet
nt        = 1000;           % plot every nt timesteps
artesian  = 1;              % allow moulin to shed water?
Qscale    = 1;              % factor to scale FOXX runoff by
E         = 10;              % enhancement factor for creep

HFdoy     = 99999999999999;%  % Prescribe an annual date of hydrofracture?  # if yes. Really high # if no.   165; % Mid June
%% set the vertical model components
dz        = 1; %  vertical spacing, meters
z         = (0:dz:H)';
time.z    = z; %save the z profile in time
%% set the duration of the model run
sec       = 86400*numofdays;   %seconds * days
dt        = 300;        % Timestep, seconds (30 minutes)
tmax      = sec;        % seconds (5 years here)
time.t    = dt:dt:tmax; % seconds
time.dt   = dt;
time.sec  = sec;
%% set Qin 
% Construct an approximate Qin for each timestep, based on air temps:
% Qin = double(Tair>C.To) .* (Tair-C.To)*const;
% Qin = zeros(size(time.t));
% Use predetermined Qins of various types
load Qcosines.mat %1 = time, 2 cosine function, 3 
Qcos2   = Qcos2(1:end,:) ;
%change Qcos2 column for different types: 1. cosine, 2. cosine with small
%melt event, 3. cosine with large melt event, 4. quasi-real data
Qin     = interp1(Qcos2(:,1), Qcos2(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
Qin     = Qin*0.5 +4; %scale Qin to deal with a few model issues
time.Qin = Qin;  %save for future plotting

% figure
% hold on
% plot(time.t(1:48*4), Qcos2(1:48*4,2))
% plot(time.t(1:48*4), Qcos2(1:48*4,3))
% plot(time.t(1:48*4), Qcos2(1:48*4,4))
% plot(time.t(1:48*4), Qcos2(1:48*4,5))
clear Qcos2
%% set Ice temperature characteristics 
Tz      = importTz('Ryser_foxx',z);
Tfar    = Tz; % Kelvin
xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx]...
        = setupx(dt,chebx,xmax);
T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1)  = C.T0;   % Melting point at the moulin wall
time.icetemp = Tz; %just save in the time file for reference
%% define initial moulin characteristics

%hw      = zeros(1,length(time.t));
hwint   = H ; %set the inital water level as 
%hw(1)   = hwint;
Mrmin   = 1e-9;  % 1 mm
M.r(:,1) = R0*ones(size(z));

%create a non cylinderical initial radius
initrad = (z+(H/0.5)) ./ (H/1);
M.r(:,1) = initrad; %To use this, the moulin should be filled 

% initalize the horizontal coordinate system
M.xu = zeros(length(z),1);
M.xd = zeros(length(z),1) + M.r.*2;
%% Set turbulence parameters

relative_roughness = 0.2; %increasing this value increases the amount of melting due to turbulence.

include_ice_temperature = true; %true means that the change in the ice temperature is included in...
%the calculated change in moulin radius. If false, it makes the implicit
%assumption that the ice temperature and water temperature are both at the pressure melting temperature. 

if include_ice_temperature
    Ti = Tz;
else
    Ti = NaN; %#ok<UNRCH>
end

Bathurst = true; %true means that the friction factor is calculated using..
%the Bathurst formulation... this equation is valid when
%roughness height ./ hydrualic diameter >= 0.05
% if false, the Colebrook-White formulation will be applied, which is only
% valid when roughness height ./ hydrualic diameter < 0.05

%% Assign elastic deformation parameters
sigx = -50e3;%100e3;
sigy = -50e3;%-100e3;
tauxy = 100e3;%100e3;


%% save general parameters in time file 
time.parameters.sigxytauxy = [sigx, sigy, tauxy];
time.parameters.relative_roughness = relative_roughness;
time.parameters.creepenhancement = E;
time.parameters.H = H;
time.parameters.L =L;
time.parameters.R0 = R0;
time.parameters.numofdays =  numofdays;
time.parameters.f = f;
%% Set up initial figure

% figure(3); clf;
% NP = 5; ii=1;
% subplot(1,NP,1); hold on; title('Moulin radius'); xlabel('meters'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Creep'); xlabel('m / dt'); ii=ii+1;
% % subplot(1,NP,ii); hold on; title('Refreezing'); xlabel('m / dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Turbulent melt'); xlabel('m /dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Potential energy melt-out'); xlabel('m /dt'); ii=ii+1;
% subplot(1,NP,ii); hold on; title('Elastic'); xlabel('m / dt')
% %
% % Record minimum and maximum moulin diameter
% time.Mrmm = zeros(2,numel(time.t));
% time.V    = zeros(size(time.t));
% time.hw   = zeros(size(time.t));
% time.Vcapacity = zeros(size(time.t));
% time.Vfrz = zeros(size(time.t));
% time.Vturb = zeros(size(time.t));
% dF = zeros(size(z));
% dE = zeros(size(z));
% dC = zeros(size(z));
% dM = zeros(size(z));
% dP = zeros(size(z));


%% Step through time
cc = 0;

for t = time.t
    
    cc = cc+1;
    % Consider using the previous moulin radius in all calculations in each
    % timestep, so that the final result is not dependent on the order in
    % which I do creep, refreeze, turbulent melt, elastic, etc.
    Mrprev = M.r;
    
    % Check if today is the day that we hydrofracture and reopen the bottom
    % of the moulin:
    if ~mod(HFdoy*86400 - t,sec)
        fprintf('Hydrofracture event! at t=%1.2f years (cc=%d)\n',t/sec,cc)
        M.r = max(M.r,R0);
    end
    %
   
%%%%%%%%%%
%Water level and subglacial conditions
%Dont let the water level go over the top of the ice    
%     if hw > H
%         time.excess(cc) = 
%         hw = H;   
%     elseif hw <= H
%         hw = hw;     
%     end

    tspan = [t,t+dt];
    y0 = [hw, S];
    %[hw,S,Qout]   = subglacialsc(Mrprev,z,Qin(cc),H,L,C,tspan,y0);
    opt = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));
    [hw,S,Qout]   = subglacialsc(Mrprev,z,Qin(cc),H,L,C,tspan,y0, opt); %consider adding Vadd to the qin values
    
    time.S(cc)    = S;
    time.hw(cc)   = hw;
    time.Qout(cc) = Qout;
    
    
% % % % %     if tmp > limit
% % % % %         break
% % % % %     end
%     % Moulin water volume:
%     V = watervolume(V,Vturb,Vfrz,Qin(cc),Qout(cc),dt);
%     % Make an artesian spring if V > Vmoulin
%     if artesian, V = min(trapz(z,pi*Mr.^2),V); end
%     % How high does that fill the moulin?
%     hw = waterlevel(Mr,z,V);
%             time.hw(cc) = hw;
%             time.V(cc) = V;



%%%%%%%%% dC: Creep deformation
%Creep deformation: do this first because it is a larger term  
    dC = creep(Mrprev,z,H,hw,T,dt,E,C);
    time.dC(:,cc) = dC;
    
    
%%%%%%%%% dF: Refreezing
% Refreezing
%     T(z>hw,1) = Tair(cc);
%     [~,dF,T,Vfrz] = refreeze(Mrprev,T,z,hw,dF,nx,x,dx,dt);
%             time.Vfrz(cc) = Vfrz;



%%%%%%%%% dM: Turbulent melting
% Turbulent melting: 
  [dM, uw, Vadd] = turbulence(hw, Qout, Mrprev, dt, Ti, z, relative_roughness, Bathurst, include_ice_temperature);
   time.dM(:,cc)  =  dM;
   time.uw(:,cc)  =  uw;
   time.Vadd(cc)  = Vadd;

%%%%%%%%%   
    %deal with the Vadd term by adding it to the next Qin timestep so that
    %it is integrated 
    if cc < length(time.t)
        Qin(cc +1) = Qin(cc+1) + Vadd./dt;
    end
    
%%%%%%%%% dE: Elastic deformation   
% Elastic deformation: do this last because it is a function of moulin 
  % radius.  Elastic deformation is small and sensitive to water pressure
    dE = elastic(z,Mrprev,hw,H,sigx,sigy,tauxy,C);
    time.dE(:,cc) = dE;

%%%%%%%%% dP: Expansion from gravitational potential energy above the water
%%%%%%%%% line
    dP = potentialdrop(Qin(cc),z,hw,Mrprev,dt,C,f);
    % The reason for calculating the above is to offset the elastic closure
    % at the top of the moulin.  On its own, the moulin will close
    % elastically after some days to months (depending on C.E).  We know
    % that does not happen.  Hence, we add some turbulent "waterfall"
    % melting above the water line.
    time.dP(:,cc) = dP;
        
    % Calculate the horizontal position of the moulin within the ice column
    M.xu = M.xu - dC - dE - dM - dP; % - dOC - dG;
    M.xd = M.xd + dC + dE + dM + dP; % +dG;
    
    % Now use the moulin positions to calculate the actual radius:
    M.r = (M.xd - M.xu) / 2;
    %M.r = M.r + dC + dE + dM + dP; % + dF + dM  + dP;
    %M.r = max(M.r,M.rmin);
        
    % Record the used moulin geometry 
    time.M.r(:,cc) = M.r;
    time.M.xu(:,cc) = M.xu;
    time.M.xd(:,cc) = M.xd;
    %
    % Record volume capacity of moulin
    time.Vcapacity(cc) = trapz(z,pi*M.r.^2);
    %
    % Plot profiles sometimes
%     if ~mod(cc,nt)
%         figure(3); ii=1;
%             [~,jw] = min(abs(z-hw));
%             subplot(1,NP,1); plot(Mr,z,Mr(jw),hw,'*k'); title(sprintf('Moulin radius at t=%1.1f days',t/3600/24)); ii=ii+1;
%             subplot(1,NP,ii); plot(dC,z,dC(jw),hw,'*k'); ii=ii+1;
% %             subplot(1,NP,ii); plot(dF,z,dF(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dM,z,dM(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dP,z,dP(jw),hw,'*k'); ii=ii+1;
%             subplot(1,NP,ii); plot(dE,z,dE(jw),hw,'*k');
%     end
    %
end
%% figures

if plot_using_lauren
    laurensplots(time, save_figures, save_location, datetime, visible_figures)
end




if save_timevariable
  cd(save_location)
  tmp = datestr(now,'mm-dd-yyyy');
  mkdir(tmp);
  cd(tmp);
  filename = ['modelrun', '_R0-', num2str(R0), '_H-', num2str(H), '_', num2str(numofdays), 'd_',  datetime, '_outputs.mat']
  save(filename, 'time')
end 


% figure(30); clf; 
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
% title(sprintf('Volume entering and leaving per dt (Q$_{in,out}$*dt), %1.2f day lag',ndaylag))
% ylabel('(m$^3$)')
% legend('Q$_{in}$ *dt','Q$_{out}$ *dt')%,'V in moulin')
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
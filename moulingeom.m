% Solve for moulin geometry, considering
% 1. elastic deformation (opening / closure)
% 2. creep deformation (opening / closure)
% 3. refreezing (closure)
% 4. turbulent melting (opening)
% 5. open-channel melting (opening)
% 6. xz-shear deformation
%
%

clear variables
close all

%Warning: Failure at t=2.700000e+04.  Unable to meet integration tolerances without reducing the step size
%below the smallest value allowed (5.820766e-11) at time t. %These warnings off is not ideal, 
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
make_simple_plots =  true;
save_figures      =  false;
visible_figures   =  true;

% Do you want to save the 'time' structure? This has (almost) everything
% documented...
save_timevariable =  false; 

save_location = './modeloutputs';
datetime = datestr(now,'mm-dd-yyyy_HHMM'); %This will assign a unique date and time for both the figures and the model outputs
%% define some basic parameters
C         = makeConstants;  %constants used for parameterizations 
Tdatatype = 'Ryser_foxx';   %ice temperature profile to extrapolate from
numofdays = 20;             %set the number of days for the model run
H         = 800;            % ice thickness, meters
R0        = 5;              % radius of moulin initially
L         = 12e3;           % Length of the subglacial channel
f         = 0.05;           % fraction of the potential energy used to open the top of the moulin (above water level)
alpha     = 0.03;           % regional surface slope (unitless), for use in Glen's Flow Law
n         = 3;              % flow law exponent (Glen's Flow Law)

%inital guesses for subglacial model
hw(1) = H;                  % moulin water level (m)
S(1)  = 1.5*R0;                 % subglacial channel cross sectional area (m^2)


chebx     = 0;              % chebx=1 is not working yet
nt        = 1000;           % plot every nt timesteps
% artesian  = 1;              % allow moulin to shed water?
% Qscale    = 1;              % factor to scale FOXX runoff by
E         = 5;             % enhancement factor for creep

% HFdoy     = 99999999999999;%  % Prescribe an annual date of hydrofracture?  # if yes. Really high # if no.   165; % Mid June
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
%melt event, 3. cosine with large melt event, 4. quasi-real data,
%5. realistic, and 6. realistic but tapering to zero without massive
%diurnal variability.

% Small melt event input:
%Qin     = interp1(Qcos2(:,1), Qcos2(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
% Big melt event input:
Qin     = interp1(Qcos2(:,1), Qcos2(:,3), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
% Quasi real/random input:
Qin     = interp1(Qcos2(:,1), Qcos2(:,5), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
Qin     = Qin*0.8 +3; %scale Qin to deal with a few model issues
time.Qin = Qin;  %save for future plotting
clear Qcos2

%% set Ice temperature characteristics 

Tfar    = importTz('Ryser_foxx',z); % Kelvin
xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx]...
        = setupx(dt,chebx,xmax,C);
T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1)  = C.T0;   % Melting point at the moulin wall
time.icetemp = Tfar; %just save in the time file for reference
%% define initial moulin characteristics

%hw      = zeros(1,length(time.t));
hwint   = H ; %set the inital water level as 
%hw(1)   = hwint;
Mrmin   = 1e-9;  % 1 mm
M.r     = R0*ones(size(z));

%create a non cylinderical initial radius
initrad   = M.r; %(z+(H/0.5)) ./ (H/1);
M.r_minor = initrad; %To use this, the moulin should be filled 
M.r_major = initrad; %To use this, the moulin should be filled 

% initalize the horizontal coordinate system
%This assumes that ice flow is from left to right 
M.xu = -M.r_major;
M.xd =  M.r_minor;
% Pin the bed of the upstream wall to x=0 while retaining the initial
% moulin shape / radius:
x0 = M.xu(1);
M.xu = M.xu - x0;
M.xd = M.xd - x0;

%% Set turbulence parameters

relative_roughness = 0.2; %increasing this value increases the amount of melting due to turbulence.
relative_roughness_OC = 1e-9;%1e-12;  % This one modifies the melt from open channel flow.

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

%And initialize the added melt components
% Vadd = 0;
% Vadd_
%% Assign elastic deformation parameters
stress.sigx = -50e3;  % compressive
stress.sigy = -50e3;  % compressive
stress.tauxy = 100e3; % shear opening

%% save general parameters in time file 
time.parameters.stress = stress;
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
    Mrminor_prev  = M.r_minor;
    Mrmajor_prev  = M.r_major;
    Mxuprev = M.xu;
    
    % which nodes are underwater or at the water line (wet) versus above the water line?
    wet = locatewater(hw,z);
   
    % Calculate hydrostatic pressures everywhere
    % Ice hydrostatic stress (INWARD: Negative)
    stress.cryo = -C.rhoi*C.g*(H-z);
    % Water hydrostatic stress (OUTWARD: Positive)
    stress.hydro = C.rhow*C.g*(hw-z);
    stress.hydro(~wet) = 0; % Anywhere that is not wet does not have the opening force from water
    
%%%%%%%%%%
% Subglacial Schoof model: Conduit size
    tspan = [t,t+dt];
    y0 = [hw, S];
            %[hw,S,Qout]   = subglacialsc(Mrminor_prev,z,Qin(cc),H,L,C,tspan,y0);
    opt = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));
            %Qin_tot       = Qin(cc) + time.V
    [hw,S,Qout]   = subglacialsc(Mrminor_prev,z,Qin(cc),H,L,C,tspan,y0, opt); %consider adding Vadd to the qin values
    
    time.S(cc)    = S;
    time.hw(cc)   = hw;
    time.Qout(cc) = Qout;
    



%%%%%%%%% dC: Creep deformation
%Creep deformation: do this first because it is a larger term  
    dC_minor = creep(Mrminor_prev,z,H,stress,T,dt,E,C);    
        time.dC_minor(:,cc) = dC_minor;    
    dC_major = creep(Mrmajor_prev,z,H,stress,T,dt,E,C);
        time.dC_major(:,cc) = dC_major;
    
    
%%%%%%%%% dF: Refreezing
% Refreezing
%     T(z>hw,1) = Tair(cc);
%     [~,dF,T,Vfrz] = refreeze(Mrminor_prev,T,z,hw,wet,dF,nx,x,dx,dt);
%             time.Vfrz(cc) = Vfrz;



%%%%%%%%% dM: Turbulent melting
% Turbulent melting: 
   [dM, uw, Vadd_turb] = turbulence(hw, Qout, Mrminor_prev,Mrmajor_prev, M.xd, dt, Ti, dz, z, wet, relative_roughness, Bathurst, include_ice_temperature);
       time.dM(:,cc)  =  dM;
       time.uw(:,cc)  =  uw;
      % time.V(cc)  = Vadd;
   %[dM_major, uw_major] = turbulence(hw, Qout, Mrmajor_prev, dt, Ti, z, relative_roughness, Bathurst, include_ice_temperature);
   %    time.dM_major(:,cc)  =  dM_major;
   %    time.uw_major(:,cc)  =  uw_major;
       %time.Vadd_major(cc)  = Vadd_major;
   % Calculate the water volume added to the moulin by calculating the enlargement of the moulin due to turbulent melting 
       %Vadd_turb = waterVolumeFromTurbulence(Mrminor_prev, Mrmajor_prev, dM, z, wet);
%        Vadd_turb = waterVolumeFromTurbulence(Mrminor_prev, Mrmajor_prev, dM, z, wet, dt);
       time.Vadd_turb(cc) = Vadd_turb;
       
%%%%%%%%% dOC: Melting due to open channel flow above the moulin water line
   [dOC, Vadd_oc] = openchannel(hw, Qin(cc), M.r_minor, M.r_major, M.xu, dt, Ti, dz, z, relative_roughness_OC, Bathurst, include_ice_temperature, wet);
  
   % Scale the open channel displacement down by 1/2 to reflect the
   % displacement at exactly the upstream point:
   dOC = dOC / 2;
   
   time.dOC(:,cc)  =  dOC;
   time.Vadd_oc(cc)    =  Vadd_oc;

%%%%%%%%% Vadd: Added water from melted ice into Qin
    % NOTE 2 MARCH 2020: This is a large amount of meltwater (~8 m2 per dt)
    % compared to the current Qin (~4 m2 per dt).  It can break the model
    % if the moulin and subglacial conduit aren't big enough.
    % Add the  
    if cc < length(time.t)
        Qin(cc+1) = Qin(cc+1) + Vadd_turb / dt + Vadd_oc / dt;
    end
    
%%%%%%%%% dE: Elastic deformation   
% Elastic deformation: This is small, and sensitive to water pressure
    dE_minor = elastic(Mrminor_prev,stress,C);
        time.dE_minor(:,cc) = dE_minor;
    dE_major = elastic(Mrmajor_prev,stress,C);
        time.dE_major(:,cc) = dE_major;

%%%%%%%%% dP: Expansion from gravitational potential energy above the water
%%%%%%%%% line
    dP = potentialdrop(Qin(cc),wet,Mrminor_prev,dt,C,f);
    % The reason for calculating the above is to offset the elastic closure
    % at the top of the moulin.  On its own, the moulin will close
    % elastically after some days to months (depending on C.E).  We know
    % that does not happen.  Hence, we add some turbulent "waterfall"
    % melting above the water line.
    time.dP(:,cc) = dP;
        
    
%%%%%%%%% dG: Asymmetric deformation due to Glen's Flow Law
    dG = deformGlen(H, T, alpha, z, n, dt, C);
    time.dG(:,cc) = dG;
    
    
    % Calculate the horizontal position of the moulin within the ice column
    M.xu = M.xu - dC_major - dE_major - dM + dG - dOC;% - 0*dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
                                % Important Note: the +dG above is correct.
                                % The upstream wall moves downstream.
                                
    M.xd = M.xd + dC_minor + dE_minor + dM + dG;% + dP;
                                % The downstream wall also moves downstream
                                % at the same rate, dG.

    % Shift them both back upstream so that the bed of the upstream wall
    % stays pinned at x = 0:
    x0 = M.xu(1);
    M.xu = M.xu - x0;
    M.xd = M.xd - x0;
    %
    % Now use the moulin positions to assign the major and minor radii:
    M.r_minor = max(M.r_minor + dC_minor + dE_minor + dM, Mrmin);
    M.r_major = (M.xd - M.xu) - M.r_minor;
        
    % Record the used moulin geometry 
    time.M.r_minor(:,cc) = M.r_minor;
    time.M.r_major(:,cc) = M.r_major;
    time.M.xu(:,cc) = M.xu;
    time.M.xd(:,cc) = M.xd;
    %
    % Record volume capacity of moulin 
    % This reflects the semi-ellipse, semi-circular geometry 
    % and uses the variable "wet" for where is water.
    [time.Mcapacity(cc), time.Wvolume(cc)] = moulincapacity(M,z,wet);

end
%% figures

if make_simple_plots
    simpleplots(time, save_figures, save_location, datetime, visible_figures)
end

%%


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
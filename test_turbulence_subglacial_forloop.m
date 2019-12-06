% Stand alone test turbulence parameters 
clear variables 
%close all
clc

C       = makeConstants;
chebx   = 0;  % chebx=1 is not working yet
dz      = 1;
H       = 1000;
z       = (0:dz:H)';

sec     = 86400*15;
dt      = 300; % seconds 300 = 5min
tmax    = sec; 
time.t  = dt:dt:tmax; 

load Qcosines.mat
Qcos2  = Qcos2(1:end,:) ;
Qin         = interp1(Qcos2(:,1), Qcos2(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
%Qin(24761:end) = 0;
Qin =Qin*0.8 +3; %for fake values
%Qin =0.8*Qin+4.5; %for real values 
figure; plot(Qin); hold on; 

%Qin(:) =3;
%The commented info below is for when the subglacial component is not
%included
% Qin(1)      = 5;
% ndaylag     = 1/24; 
% nlag = round(ndaylag*24*3600/dt);
% Qout = fastsmooth(Qin,nlag/2,3,1);
% % Qout = max(0,Qout);
% Qout = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
% Qout = Qout * 1;

%  Elastic deformation parameters %increase to reduce elastic response
sigx = -50e3;%100e3;
sigy = -50e3;%-100e3;
tauxy = 100e3;%100e3;

% Assign hydraulic gradient parameters and water flow speed:
nz = numel(z);
ubottom = 1;  % m/s; this is something like 1 mm/sec
utop = ubottom; % just make something up for now
u = linspace(ubottom,utop,nz)';
u0 = ubottom; z0 = 0;
L = 15e3; % length scale over which to take hydraulic gradient

%% create a non cylinderical initial radius
 initrad = (z+(H/0.5)) ./ (H/1);
 
%%

ks = 0.1;

include_ice_temperature = true; %true means that the change in the ice temperature is included in...
%the calculated change in moulin radius. If false, it makes the implicit
%assumption that the ice temperature and water temperature are both at the pressure melting temperature. 

Bathurst = true; %true means that the friction factor is calculated using..
%the Bathurst formulation... this equation is valid when
%roughness height ./ hydrualic diameter >= 0.05
% if false, the Colebrook-White formulation will be applied, which is only
% valid when roughness height ./ hydrualic diameter < 0.05

% create ice temperature profiles
if include_ice_temperature
    Ti = importTz('Ryser_foxx', z);
    %Ti      = flipud(Ti);
else
    Ti = NaN; %#ok<UNRCH>
end

Tz      = importTz('Ryser_foxx',z);
%Tz      = flipud(Tz);
Tfar    = Tz; % Kelvin
xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx]...
        = setupx(dt,chebx,xmax);
T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1)  = C.T0;   % Melting point at the moulin wall
E       = 10; %enhancement factor


%initialize variables 

hw      = zeros(1,length(time.t));
hwint   = H;
hw(1)   = hwint;
R0      = 2.5;  % radius of moulin initially
%Mr(:,1) = R0*ones(size(z));
Mr(:,1) = initrad; %To use this, the moulin should be filled 
MoulinSysVol(:,1) ...
        = pi .* Mr(:,1).^2;

S       = nan .* zeros(1, length(time.t));
S(1)    = 2;
% Qout    = nan .* zeros(1, length(time.t));
% Qout(1) = Qin(1);
SubSysVol(1) ...
        = S(1) .* L;

 %TotWaterVol(1) ...
 %        = SubSysVol(1) + hw(1).*pi .* (R0.^2);
TotWaterVol(1) ...
        = SubSysVol(1) + hw(1).*pi .* (mean(Mr(:,1))).^2;

ExVol(1)= 0;    

% do a for loop to calculate turbulence

for ii = 1:length(time.t)
    
    % calculate the Qout for a given subglacial configuration and moulin
    % head 
    Qout(ii)  = C.c3 .* (S(ii).^(5/4)) ...
               .* (((C.rhow .* C.g .* hw(ii)) ./L).^ (0.5)); %There may be times i.e. when hw > 0.91*H, where there should be some type of sheet flow

           
           %     tmp = sqrt((S(ii) ./pi));
%     psi = C.rhow .* C.g .*hw(ii) ./L;
%     Qout1(ii) = S(ii) .* (((((pi+2).* psi .*tmp) + (2.* tmp .* psi))./(C.rhow .* C.f)).^(0.5));
      %Qout(ii) = sqrt(((S(ii).^(8/3)) .* psi)./ 850);
            if Qout(ii) < 0 %trigger only if something odd has happened and Qout becomes negative
                Qout(ii)        = 0;
                 disp('Qout is negative!!')
                 disp('Assigning Qout==0')
            end
            
    waterpresent = hw(1,ii) - z; %logical matrix to determine in what nodes water is present
    waterpresent(waterpresent<0) =0;
    
    S_moulin(:,ii) = (pi * Mr(:,ii).^2);
        
    uw_tmp  = Qout(ii) ./  S_moulin(:,ii); %calculate the water velocity in each node
    uw_tmp(waterpresent ==0) =0; % if there is no water in a given cell,
    
    % % ------------------------------
    % % Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity.
    % %It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross
    % % section.
    if uw_tmp > 9.3
        disp('big velocity !!!')
        disp('assigning terminal velocity of 9.3ms-1')
    end
    
    uw_tmp = min(uw_tmp,9.3);
    % % ------------------------------
    uw(:,ii) = uw_tmp;
    
    Mp(:,ii)   = 2 .* pi .* Mr(:,ii); % wetted/melting perimeter
    Dh(:,ii)   = (4*(pi .* Mr(:,ii).^2)) ./ Mp(:,ii); %hydrualic diameter
    Rh(:,ii)   = (pi.* Mr(:,ii).^2) ./ Mp(:,ii); % hydraulic radius
    
    
    % calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
    Pw(:,ii)   = C.rhow .* C.g .* waterpresent;
    Tmw(:,ii)  =  273.16 - 9.8e-8.*(Pw(:,ii) - 611.73);
    
    % select the appropriate parameterization of DW friction factor
    if Bathurst
        fR(:,ii) =10* 1./((-1.987.* log10(ks./(5.15.*Rh(:,ii)))).^2); %#ok<UNRCH> %Bathurst parameterization for DW friction factor
    else
        fR(:,ii) = (1./(-2.*log10((ks./Dh(:,ii))./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
        %fR(:,ii) = 0.01;
    end
    
    dz          = nanmean(diff(z));
    % calculate head loss following Munson 2005
    hL(:,ii)  =  ((uw(:,ii).^2) .* fR(:,ii) .* dz) ./(2 .* Dh(:,ii) .* C.g);
    
    if include_ice_temperature
        dM_dt(:,ii) =    (C.rhow .* C.g .* Qout(ii) .* (hL(:,ii)./dz)) ...
            ./ (2 .* pi .* Mr(:,ii) .* C.rhoi .* (C.cw .* (Tmw(:,ii)- Ti ) + C.Lf));
        %This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
        % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
        % include the surrounding ice temperature
        
    else
        dM_dt(:,ii) =   (C.rhow .* C.g .* Qout(ii) .* (hL(:,ii)./dz)) ./ (2 .* pi .* Mr(:,ii) .* C.rhoi .* C.Lf); %#ok<UNRCH>
        % This parameterization is closer to that which is traditionally used
        % to calculate melting within a subglacial channel where the ice and
        % water are at the same temperature
    end
    
    % Make sure that there is no melting in places without water
    dM_dt_tmp = dM_dt(:,ii);
    dM_dt_tmp(waterpresent == 0) = 0;
    dM_dt(:,ii) = dM_dt_tmp;
    
    %This is where something needs to be added for melting above the water
    %line...
    
    
    dM(:,ii)  = dM_dt(:,ii) .* dt; %change in radius over the given time step
    
    Vadd(ii) = C.rhoi/C.rhow * trapz(2*pi*Mr(:,ii).*dM(:,ii), z);
    
    dC(:,ii) = creep(Mr(:,ii),z,H,hw(ii),T,dt,E,C);
    dE(:,ii) = elastic(z,Mr(:,ii),hw(ii),H,sigx,sigy,tauxy,C);
%     if 
%         dE(:,ii) =0;
%         disp('dE = 0 triggered')
%     end
    % dealing with some elastic issues
    dE_tmp = dE(:,ii);
    dE_tmp(waterpresent == 0) = 0;
    %dE_tmp(dE_tmp>0) = 0;    
    dE(:,ii) = dE_tmp;
   
    
    %subglacial system
    dS(ii)    = dt.* (C.c1 .* Qout(ii) .* ((C.rhow .* C.g .* hw(ii)) ./L)...
                 - C.c2 .* ((C.rhoi .* C.g .* H - C.rhow .* C.g .* hw(ii)).^C.n) .* S(ii)); %change in subglacial channel geometry
        
    % advance the for loop
    %1/calculate the new moulin radius
    Mr(:,ii+1) = Mr(:,ii) + dC(:,ii)+ dM(:,ii) + dE(:,ii); % moulin radius
    %2/calculate the new subglacial cross-sectional area
    S(ii+1)    = S(ii) + dS(ii);
    
    %3/calculate the volume of the new system
    SubSysVol(ii+1)  =  S(ii+1).*L; %thenew cross sectional area times the length
    MoulinSysVol(:,ii+1) = pi .* (Mr(:,ii+1).^2) .* dz; %new radius *pi *node length (depth and time field)
    
    %4/calculate new volume of water in the system
    TotWaterVol(ii+1) = TotWaterVol(ii) + Qin(ii).* dt - Qout(ii).*dt + Vadd(ii) + dS(ii).*(C.rhoi./C.rhow);
    %total vol for next timestep
                  %prev Total vol + Qin over time step - qout over time
                  %step + volume of water created from turbulent melting in the moulin (from melting ice)
                  % + volume of water created from turbulent melting in
                  % subglacial channel (from melting ice)
    
    %5/calculate new moulin water level 
    watertofillmoulin = TotWaterVol(ii+1) - SubSysVol(ii+1); %calculate the volume of water remaining to fill the moulin
%     if watertofillmoulin<0
%         watertofillmoulin =0;
%     end
    
    mv_tmp            = cumsum(MoulinSysVol(:,ii+1)); %calculate the cumulative volumen of water
    [m0, min_index]   = min(abs(watertofillmoulin-mv_tmp)); %find the nearest node/height
    extra_water       = watertofillmoulin-mv_tmp(min_index); %determine if a small up or down adjustment needs to be made
    
    
    hw(ii+1)   = z(min_index) + extra_water ./(pi .* (Mr(min_index,ii+1).^2) ); %calculate the small adjustment
    
    if hw(ii+1) > H
        hw(ii+1)    = H;
        ExVol(ii+1) = TotWaterVol(ii+1) - SubSysVol(ii+1) - sum(MoulinSysVol(:,ii+1));
        %disp('Water above ice surface !!!')
        
    end
    
    if hw(ii+1) < 0
        hw(ii+1) =0;
        disp('Something bad happened: hw is less than zero !!!')
        disp('Perhaps a timestep issue')
    end

end

%%

figure
subplot(2,1,1)
hold on
plot([0 3000], [H*0.91 H*0.91])
plot(hw)
title('hw, straight line is flotation')
subplot(2,1,2)
hold on
plot(Qin)
plot(Qout)
legend('Qin', 'Qout')
%%
spacing = 100;
endlength = 4320-1;
color1 = brewermap(endlength, 'spectral');

figure
subplot(1,4,1)
hold on
for ii = 1:spacing:endlength
plot(Mr(:,ii), z, 'color', color1(ii,:))
end
%plot([0 3000], [700*0.91 700*0.91])
%axis([0 3 0 H])
title('Radius')

subplot(1,4,2)
hold on
for ii = 1:spacing:endlength
plot(dM(:,ii), z, 'color', color1(ii,:))
end
title('Melting')



subplot(1,4,3)
hold on
for ii = 1:spacing:endlength
plot(dC(:,ii), z, 'color', color1(ii,:))
end
title('Creep')

subplot(1,4,4)
hold on
for ii = 1:spacing:endlength
plot(dE(:,ii), z, 'color', color1(ii,:))
end
title('Elastic')



%%
figure; 
set(gcf, 'position', [ 1006         349        2093         981])

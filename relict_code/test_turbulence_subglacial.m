% Stand alone test turbulence parameters 
clear variables 
close all
clc

C       = makeConstants;
chebx   = 0;  % chebx=1 is not working yet
dz      = 1;
H       = 700;
z       = (0:dz:H)';

R0      = 2;  % radius of moulin initially
Mr      = R0*ones(size(z));

sec     = 86400*10;
dt      = 900; % 1 minutes
tmax    = 0.5*sec; 
time.t  = dt:dt:tmax; 

load Qsine.mat
Qsine  = Qsine(8:end,:) ;
Qin         = interp1(Qsine(:,1), Qsine(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
Qin(Qin>5) =5;
Qin =Qin*0.5 +1.5;
%The commented info below is for when the subglacial component is not
%included
% Qin(1)      = 5;
% ndaylag     = 1/24; 
% nlag = round(ndaylag*24*3600/dt);
% Qout = fastsmooth(Qin,nlag/2,3,1);
% % Qout = max(0,Qout);
% Qout = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
% Qout = Qout * 1;

%  Elastic deformation parameters
sigx = -50e3;%100e3;
sigy = -50e3;%-100e3;
tauxy = 100e3;%100e3;

% Assign hydraulic gradient parameters and water flow speed:
nz = numel(z);
ubottom = 1;  % m/s; this is something like 1 mm/sec
utop = ubottom; % just make something up for now
u = linspace(ubottom,utop,nz)';
u0 = ubottom; z0 = 0;
L = 20e3; % length scale over which to take hydraulic gradient

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
    Ti = importTz('HarrS4C', z);
else
    Ti = NaN; %#ok<UNRCH>
end

Tz      = importTz('HarrS4C',z);
Tfar    = Tz; % Kelvin
xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
[x,dx,nx]...
        = setupx(dt,chebx,xmax);
T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
T(:,1)  = C.T0;   % Melting point at the moulin wall
E       = 1; %enhancement factor


%initialize variables 

hw      = zeros(1,length(time.t));
hwint   = H*.91-100;
hw(1)   = hwint;

MVol    = zeros(1,length(time.t));
MVol(1) = hw(1) .* (pi .* Mr(1,1).^2) ;

S       = nan .* zeros(1, length(time.t));
S(1)    = 2;
Qout    = nan .* zeros(1, length(time.t));
Qout(1) = Qin(1);



% do a for loop to calculate turbulence

for ii = 1:length(time.t)-1
    
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
    
    
    dM(:,ii)  = dM_dt(:,ii) .* dt; %change in radius over the given time step
    
    Vadd(ii) = C.rhoi/C.rhow * trapz(2*pi*Mr(:,ii).*dM(:,ii), z);
    
    dC(:,ii) = creep(Mr(:,ii),z,H,hw(ii),T,dt,E,C);
    dE(:,ii) = elastic(z,Mr(:,ii),hw(ii),H,sigx,sigy,tauxy,C);
    
    
    
    
    %advance the forloop
    Mr(:,ii+1) = Mr(:,ii) + dM(:,ii) + dC(:,ii) + dE(:,ii); % moulin radius
    MVol(ii+1) = MVol(ii) +  Qin(ii) .*dt - Qout(ii).*dt + Vadd(ii); %moulin volume
    
    % calculate Qout using Celia's function
    % Find the water level in the moulin at this timestep
   % tspan = [time.t(ii),time.t(ii) + dt];
        tspan = [0, dt];

    %     if ii == 1 %this provides the initial guess for the moulin water level
    %                 % and Qout for the subglacial component
    %         y0=[hwint 0.5];
    %     else
    %         y0 = [hw(ii), S(ii)];
    %     end
    
    y0 = [hw(ii), S(ii)];
    
    [hw(ii+1), S(ii+1), Qout(ii+1)] = subglacialsc(Mr(:,ii+1),z,Qin(ii+1),H,L,C,tspan,y0);
    
    if hw(ii+1) > H
        hw(ii+1)    = H;
        
    elseif isnan(hw(ii+1))
        hw(ii+1)    = H;
        Qout(ii+1)  = C.c3 .* (S(ii).^(5/4)) .* (((C.rhow .* C.g .* H) ./L).^ (0.5)); %some estimate of sheet opening (?)
        S(ii+1)     = S(ii) + C.c1 .* Qout(ii) .* ((C.rhow .* C.g .* H) ./L)...
            - C.c2 .* ((C.rhoi .* C.g .* H - C.rhow .* C.g .* H).^C.n) .* S(ii); % Schoof SI equation 1 without sliding opening
    else
        hw(ii+1)    = hw(ii+1);
    end
    
    if Qout(ii+1) < 0
        Qout(ii+1) = 0;
        ii
    end
    
    % for determing new water level when
    % mv_tmp = pi .* (Mr(:,ii+1).^2) .* dz;
    % mv_tmp = cumsum(mv_tmp);
    % [m0, min_index] = min(abs(MVol(ii+1)-mv_tmp));
    % extra_water =  MVol(ii+1)-mv_tmp(min_index);
    %
    %
    % hw(ii+1)   = z(min_index) + extra_water ./(pi .* (Mr(min_index,ii+1).^2) );
    % hw_tmp = hw(ii+1);
    % hw_tmp(hw_tmp <0) =0;
    % hw_tmp(hw_tmp >H*.91) =H*.91;
    % hw(ii+1) = hw_tmp;
end

%%

figure
subplot(2,1,1)
hold on
plot([0 3000], [700 700])
plot(hw)
subplot(2,1,2)
hold on
plot(Qin)
plot(Qout)

%%
spacing = 1;
endlength = length(Pw);
color1 = brewermap(endlength, 'spectral');

figure
subplot(1,4,1)
hold on
for ii = 1:spacing:endlength
plot(Mr(:,ii), z, 'color', color1(ii,:))
end
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

% 
% waterpresent = hw - z; %logical matrix to determine in what nodes water is present 
% waterpresent(waterpresent<0) =0;
% S = (pi * Mr.^2);
% uw = Qout ./  S; %calculate the water velocity in each node 
% uw(waterpresent ==0) =0; % if there is no water in a given cell,
% 
% % ------------------------------
% % Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity.
% %It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross
% % section.
% if uw > 9.3
%     disp('big velocity !!!')
%     disp('assigning terminal velocity of 9.3ms-1')
% end
% 
% uw = min(uw,9.3);
% % ------------------------------
% 
% Mp   = 2 .* pi .* Mr; % wetted/melting perimeter
% Dh   = (4*(pi .* Mr.^2)) ./ Mp; %hydrualic diameter
% Rh   = (pi.* Mr.^2) ./ Mp; % hydraulic radius
% % just give a flag if the ks/dh ratio is less than 0.05
% if (ks/Dh) < 0.05
%     disp('ks/dh < 0.05')
%     disp('Should be using Colebrook-White')
% end
% 
% 
% 
% 
% % select the appropriate parameterization of DW friction factor
% if Bathurst
%     fR = 1./((-1.987.* log10(ks./(5.15.*Rh))).^2); %Bathurst parameterization for DW friction factor
% else
%     fR = (1./(-2.*log10((ks/Dh)./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
% end
% 
% 
% % calculate the lengthscale over which head loss is determined...
% % this is simply the length of the model grid cells unless they become
% % non-uniform
% dz          = nanmean(diff(z));
% 
% % calculate head loss following Munson 2005
% hL  =  ((uw.^2) .* fR .* dz) ./(2 .* Dh .* C.g);
% 
% 
% if include_ice_temperature
%     dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dz)) ./ (2 .* pi .* Mr .* C.rhoi .* (C.cw .* (Tmw - Ti) + C.Lf)); 
%     %This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013)
%     % Gulley et al. (2014), Nye (1976) & Fountain & Walder (1998) to
%     % include the surrounding ice temperature 
%     
% else
%     dM_dt =    (C.rhow .* C.g .* Qout .* (hL./dz)) ./ (2 .* pi .* Mr .* C.rhoi .* C.Lf); %#ok<UNRCH>
%     % This parameterization is closer to that which is traditionally used
%     % to calculate melting within a subglacial channel where the ice and
%     % water are at the same temperature
% end
% 
% % Make sure that there is no melting in places without water
% dM_dt(waterpresent == 0) = 0;
% 
% dM  = dM_dt .* dt; %change in radius over the given time step
% 
% Vadd = C.rhoi/C.rhow * trapz(2*pi*Mr.*dM, z); %volume of meltwater gained due to melting the surrounding ice
% 
% 
% 
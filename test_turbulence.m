% Stand alone test turbulence parameters 
clear variables 
close all
clc

C = makeConstants;

dz = 1;
H  = 700;
z = (0:dz:H)';

R0 = 2;  % radius of moulin initially
Mr = R0*ones(size(z));

sec = 86400*120;
dt = 3600*24 * (0.0208333333333333/2); % 0.25 h or 15 minutes
tmax = 0.5*sec; 
time.t = dt:dt:tmax; 

load Qsine.mat
Qsine(:,2) = Qsine(:,2) +5;
Qin = interp1(Qsine(:,1), Qsine(:,2), time.t, 'linear', 'extrap'); % run an interp just in case the timeframe changes
Qin(1) = 5;
Qin = Qin;
ndaylag = 1/24; 
nlag = round(ndaylag*24*3600/dt);
Qout = fastsmooth(Qin,nlag/2,3,1);
% Qout = max(0,Qout);
Qout = [Qout(end-nlag+1:end) Qout(1:end-nlag)];
Qout = Qout * 1;

%%

ks = 1;

include_ice_temperature = true; %true means that the change in the ice temperature is included in...
%the calculated change in moulin radius. If false, it makes the implicit
%assumption that the ice temperature and water temperature are both at the pressure melting temperature. 

Bathurst = false; %true means that the friction factor is calculated using..
%the Bathurst formulation... this equation is valid when
%roughness height ./ hydrualic diameter >= 0.05
% if false, the Colebrook-White formulation will be applied, which is only
% valid when roughness height ./ hydrualic diameter < 0.05


if include_ice_temperature
    Ti = importTz('HarrS4C', z);
else
    Ti = NaN; %#ok<UNRCH>
end

%initialize variables 

hw = zeros(1,length(time.t));
hw(1) = H*.91;

MVol = zeros(1,length(time.t));
MVol(1) = hw(1) .* (pi .* Mr(1,1).^2) ;
% do a for loop to calculate turbulence

for ii = 1:length(time.t)
    
waterpresent = hw(1,ii) - z; %logical matrix to determine in what nodes water is present 
waterpresent(waterpresent<0) =0; 
    
S(:,ii) = (pi * Mr(:,ii).^2);


uw_tmp  = Qout(ii) ./  S(:,ii); %calculate the water velocity in each node 
uw_tmp(waterpresent ==0) =0; % if there is no water in a given cell, 

% % ------------------------------
% % Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity.
% %It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross
% % section.
% if uw_tmp > 9.3
%     disp('big velocity !!!')
%     disp('assigning terminal velocity of 9.3ms-1')
% end
% 
% uw_tmp = min(uw_tmp,9.3);
% % ------------------------------ 

uw(:,ii) = uw_tmp;



Mp(:,ii)   = 2 .* pi .* Mr(:,ii); % wetted/melting perimeter
Dh(:,ii)   = (4*(pi .* Mr(:,ii).^2)) ./ Mp(:,ii); %hydrualic diameter
Rh(:,ii)   = (pi.* Mr(:,ii).^2) ./ Mp(:,ii); % hydraulic radius
% just give a flag if the ks/dh ratio is less than 0.05

% calculate the pressure melting temperature of ice/water %https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
Pw(:,ii)   = C.rhow .* C.g .* waterpresent;
Tmw(:,ii)  =  273.16 - 9.8e-8.*(Pw(:,ii) - 611.73); 

% select the appropriate parameterization of DW friction factor
if Bathurst
    fR(:,ii) = 1./((-1.987.* log10(ks./(5.15.*Rh(:,ii)))).^2); %Bathurst parameterization for DW friction factor
else
    %fR(:,ii) = (1./(-2.*log10((ks./Dh(:,ii))./3.7))).^2; %#ok<UNRCH> %Colebrook-White parameterization for DW friction factor
    fR(:,ii) = 50;
end

dz          = nanmean(diff(z));
% calculate head loss following Munson 2005
hL(:,ii)  =  ((uw(:,ii).^2) .* fR(:,ii) .* dz) ./(2 .* Dh(:,ii) .* C.g);

if include_ice_temperature
    dM_dt(:,ii) =    (C.rhow .* C.g .* Qout(ii) .* (hL(:,ii)./dz)) ...
               ./ (2 .* pi .* Mr(:,ii) .* C.rhoi .* (C.cw .* (Tmw(:,ii) - Ti) + C.Lf)); 
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





%advance the forloop
Mr(:,ii+1) = Mr(:,ii) + dM(:,ii); % moulin radius
MVol(ii+1) = MVol(ii) +  Qin(ii) .*dt - Qout(ii).*dt + Vadd(ii); %moulin volume

mv_tmp = pi .* (Mr(:,ii+1).^2) .* dz;
mv_tmp = cumsum(mv_tmp);
[m0, min_index] = min(abs(MVol(ii+1)-mv_tmp));
extra_water =  MVol(ii+1)-mv_tmp(min_index);


hw(ii+1)   = z(min_index) + extra_water ./(pi .* (Mr(min_index,ii+1).^2) );
hw_tmp = hw(ii+1);
hw_tmp(hw_tmp <0) =0;
hw_tmp(hw_tmp >H) =H;
hw(ii+1) = hw_tmp;
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

figure
subplot(1,2,1)
hold on
plot(Mr, z)

subplot(1,2,2)
plot(Ti-273.15,z)
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

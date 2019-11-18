% Homework problem
%
% What are the creep closure rates in cylindrical boreholes to various
% depths,
% at various water pressures less than or equal to flotation?
%
% Compare to Naruse et al. 1988, Annals of Glaciology Vol. 11 page 100
% "Closure rate of a 700m deep borehole at Mizuho Station, East Antarctica"
%
%
% Based on DynamicDomain_CreepClosure_radial.m
%      and DynamicDomain_CreepClosure_radial_hydrostatic.m
%
clear all
G = makeGlobalParams;

pcs = zeros(11,3);
    pcs(1,:) = [0.1 0.4 1];          %blue
    pcs(2,:) = [0.7 0.15 0.15];     % marroon %[0.4 0.4 0.4];      %dark gray
    pcs(3,:) = [0.1 0.8 0.3];      %green
    pcs(4,:) = [0.1 0.9 0.9];      %cyan
    pcs(5,:) = [1 0.7 0.4];        %orange  pale orange [1 0.95 0.7]
    pcs(6,:) = [0.7 0.1 0.7];      %magenta
    pcs(7,:) = [0.8 .8 0.8];       %gray
    pcs(8,:) = [0.95 0.25 0.25];  %red
    pcs(9,:) = [0.5 .2 0.9];       %purple
    pcs(10,:) = [1 0.95 0.3];      %yellow/gold
    pcs(11,:) = [0.95 0.55 0.9];    %pink
%
pressure = 0;%5880000;         % water pressure in the borehole (Pa)
Z = 700;                 % borehole height (m)
A = 2.5e-26;              % Flow law parameter (Pa^-3 s^-1, for -35*C ice)
%A = 4.4e-25;             % Flow law parameter (Pa^-3 s^-1, for -10*C ice)
alpha = 0.001;             % Surface slope of ice
sigmaT = 100e3;           % background tensile (-) / compressive (+) stress (Pa)
D0 = 0.175 / 2;                 % Original borehole radius (m)
%
% Dates of borehole measurements (Figure 2):
%
timevec = [datenum('03-Aug-1984') datenum('23-Feb-1985') datenum('30-Apr-1985') datenum('06-Jun-1985') datenum('01-Jul-1985') datenum('02-Aug-1985') ...
           datenum('12-Sep-1985') datenum('17-Jan-1986') datenum('09-May-1986') datenum('07-Jul-1986')] * 3600*24; 
%
%
% Water pressure determines water depth:
w = Z - pressure / G.rhow / G.g;
% Vertical coordinate
z = 0:1:Z;
%
% Calculate vertical stresses
sigzi = G.rhoi*G.g*z;
sigzw = -G.rhow*G.g*max(z-w,0);
sigmaZ = sigzi + sigzw;
sigxz = sigzi * sin(alpha);  % assume surface slope alpha = 0.01
%
epsdotA = A/3 * (sigmaZ / 3).^3;
epsdotB = A/2 * sigmaT * (sigmaT.^2 + 1/4 * sigmaZ.^2 + sigxz.^2);
%
%
DA(1,:) = D0 * ones(size(z));
DB = DA;
%
% Plot hydrostatic pressure
        figure(6); clf; subplot(1,3,1,'replace'); hold on;
        plot(sigmaZ/1e6,z,'k','linewidth',2)%,'color',[0.1 0.8 0.3])
        set(gca,'ydir','reverse')
        xlabel('$\sigma_z$ (MPa)')
        ylabel('Depth in ice column (m)')
        title('Vertical stresses');
        text(min(sigmaZ/1e6),mean(z),sprintf('%2.0f MPa water pressure ',pressure/1e6),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
%
% Plot borehole diameter
        subplot(1,3,2:3,'replace'); hold on
        plot(2*DA*1e3,z,'-k')%'color',pcs(1,:));
        set(gca,'ydir','reverse')
        set(gca,'xlim',[65 175])
        xlabel('Borehole diameter (mm)')
        title('1984 Borehole at Mizuho Station, EAIS')        
        text(165,mean(z),sprintf('Strain rate: using hydrostatic stress only'),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
%{        
% Plot borehole diameter
        subplot(1,3,3,'replace'); hold on
        plot(2*DB*1e3,z,'-k')%'color',pcs(1,:));
        set(gca,'ydir','reverse')
        set(gca,'xlim',[65 175])
        xlabel('Borehole diameter (mm)')      
        text(165,mean(z),sprintf('Strain rate: using tensile stress %2.0f kPa, alpha=%0.3f',sigmaT/1e3,alpha),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
%}      
%    
%
i = 1;
t1 = timevec(1);

%
for ti = 2:length(timevec)
    
    t0 = t1;
    t1 = timevec(ti);
    dt = t1 - t0;    
    %    
    %
    % Calculate closure in this timestep
    DA(ti,:) = DA(ti-1,:) .* exp(-epsdotA * dt);
    %DB(ti,:) = DB(ti-1,:) .* exp(-epsdotB * dt);
    %
    % Plot change in crevasse wall
    subplot(1,3,2:3)
    plot(2*DA(ti,:)*1e3,z,'color',pcs(ti,:));
    %
    %subplot(1,3,3)
    %plot(2*DB(ti,:)*1e3,z,'color',pcs(ti,:));
    %
    %   
    %
end
%
%{
dt = 1e0;
DAold = D0;
for t = timevec(1):dt:timevec(end)
    DA = DAold .* exp(-epsdotA * dt);
    DAold = DA;
    if min(abs(timevec-t)) <= dt/2
        [~,ti] = min(abs(timevec-t));
        subplot(1,3,2); plot(2*DA*1e3,z,'color',pcs(ti,:));
        drawnow
        %fprintf('ti=%d\n',ti)
    end
    
end
%}
%{
dt = 1e5;
nsteps = round( (timevec - timevec(1)) / dt);
for ti = 1:length(timevec)
    D(:,ti) = D0 * exp(-epsdotA * dt *nsteps(ti));
    subplot(1,3,2:3); plot(2*D(:,ti)*1e3,z,'color',pcs(ti,:));
end
%}

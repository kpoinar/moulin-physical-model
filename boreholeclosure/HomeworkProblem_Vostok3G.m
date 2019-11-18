
% Borehole 3G at Vostok
% by Blinov and Dmitriev (1987) and Salamatin et al (1998)
% From Table 4 in Talalay and Hooke, 2007 (Annals)
% 
clear all;
G = makeGlobalParams;
pcs = zeros(11,3);
    pcs(1,:) = [0.1 0.4 1];          %blue
    pcs(2,:) = [0.7 0.15 0.15];     % marroon 
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
% The hole was filled to within 135 m of the surface with low-temperature
% aircraft fuel, TS-1, mixed with 9% (by weight) of densifier
% (trichlorofluoromethane - CFC 11) to minimize closure.
% TS-1 density: 787 kg/m3
% CFC-11 density: 1.5432 g.cm?3
% Total density: 855 kg/m3
rhofluid = 855;
%
% Depth Temperature Pressure difference (MPa) Age of ice (ka) Borehole diameter (mm)
% 19 Jul 1986
% 27 May 1988
% 31 Aug 1988
% 24 Nov 1988
% 27 Jun 1989
% 08 Oct 1989
% 04 Jan 1990
% 18 Oct 1990
data(:,1) = [1000 -49.4 -1.119 68.1 151 148 148 148 147 147 146 145];
data(:,2) = [1100 -48.4 -1.141 76.0 148 147 147 146 146 146 146 144];
data(:,3) = [1200 -47.3 -1.163 83.8 148 147 147 146 146 146 146 144];
data(:,4) = [1300 -46.1 -1.185 91.3 149 147 147 146 145 145 145 144];
data(:,5) = [1400 -44.9 -1.208 98.8 148 145 145 145 145 144 144 144];
data(:,6) = [1500 -43.6 -1.230 106.1 116 114 113 113 111 111 110 110];
data(:,7) = [1600 -42.3 -1.253 113.3 116 114 113 113 112 111 111 111];
data(:,8) = [1700 -40.9 -1.276 120.6 116 114 114 113 111 111 110 110];

dates(1) = datenum('19 Jul 1986');
dates(2) = datenum('27 May 1988');
dates(3) = datenum('31 Aug 1988');
dates(4) = datenum('24 Nov 1988');
dates(5) = datenum('27 Jun 1989');
dates(6) = datenum('08 Oct 1989');
dates(7) = datenum('04 Jan 1990');
dates(8) = datenum('18 Oct 1990');
dates = dates - dates(1);         % number of days since 19 Jul 1986
timevec = dates * 24*3600;        % seconds since then


T = mean(data(2,:));
A = G.A0less * exp(-G.Qless / G.R / (G.To + T));
z = (0:10:1700)';
w = 135;
sigzi = G.rhoi*G.g*z;
sigzw = -rhofluid*G.g*max(z-w,0);
pressure = max(0,(z-w) * rhofluid * G.g);
sigmaZ = sigzi + sigzw;
epsdot = A/5 * (sigmaZ/3).^3;
D0 = 0.15/2;
D = D0 * ones(size(z,1),size(timevec,2));
[junk,kink] = min(abs(z-1400));
D(kink:end,:) = 0.116/2;
D0 = D(:,1);
%
t1 = timevec(1);
figure(1); subplot(1,3,1,'replace')
        plot(sigmaZ/1e6,z,'k','linewidth',2)%,'color',[0.1 0.8 0.3])
        set(gca,'ydir','reverse')
        xlabel('$\sigma_z$ (MPa)')
        ylabel('Depth in ice column (m)')
        title('Vertical stresses');
        text(min(sigmaZ/1e6),mean(z),sprintf('%Borehole filled with fluid to depth of %2.0f meters ',w),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
        
subplot(1,3,2:3,'replace'); hold on
        plot(2*D0*1e3,z,'-k')%'color',pcs(1,:));
        set(gca,'ydir','reverse','yticklabel',[])
        %set(gca,'xlim',[165 175])
        xlabel('Borehole diameter (mm)')
        title('3G Borehole at Vostok')        
        text(165,mean(z),sprintf('Strain rate: using hydrostatic stress only Eps = -A/2*(1/sqrt(3)*sigZ)^3'),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
%{
for ti = 2:length(timevec)
    %
    t0 = t1;
    t1 = timevec(ti);
    dt = t1 - t0;    
    %
    % Calculate closure in this timestep
    D(ti,:) = DA(ti-1,:) .* exp(-epsdot * dt);
    %
    % Plot change in crevasse wall
    subplot(1,3,2:3)
    plot(2*D(ti,:)*1e3,z,'color',pcs(ti,:));
    %
    %
end
%}
dt = 1e-5;
nsteps = round( (timevec - timevec(1)) / dt);
for ti = 1:length(timevec)
    D(:,ti) = D0 .* exp(-epsdot * dt *nsteps(ti));
    subplot(1,3,2:3); plot(2*D(:,ti)*1e3,z,'color',pcs(ti,:));
    plot(data(ti+4,:),data(1,:),'.','color',pcs(ti,:),'markersize',40);
end
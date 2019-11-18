
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

dates = [0 79 179];         % number of days since 19 Jul 1986
timevec = dates * 24*3600;        % seconds since then


T = -0;
A = G.A0less * exp(-G.Qless / G.R / (G.To + T));
z = (0:1:150)';
w = 150;
sigzi = G.rhoi*G.g*z;
sigzw = -G.rhow*G.g*max(z-w,0);
pressure = max(0,(z-w) * G.rhow * G.g);
sigmaZ = sigzi + sigzw;
epsdot = A/2 * (1/3 * sigmaZ).^3;
D0 = 0.100/2;
D = D0 * ones(size(z,1),size(timevec,2));
D0 = D(:,1);
%
t1 = timevec(1);
figure(1); subplot(1,3,1,'replace')
        plot(sigmaZ/1e6,z,'k','linewidth',2)%,'color',[0.1 0.8 0.3])
        set(gca,'ydir','reverse')
        set(gca,'ytick',0:10:150)
        xlabel('\sigma_z (MPa)')
        ylabel('Depth in ice column (m)')
        title('Vertical stresses');
        text(min(sigmaZ/1e6),mean(z),sprintf('%Borehole filled with fluid to depth of %2.0f meters ',w),...
            'fontsize',14,'fontname','courier','fontweight','bold','color',[0 0 0])
        
subplot(1,3,2:3,'replace'); hold on
        plot(2*D0*1e3,z,'-k')%'color',pcs(1,:));
        set(gca,'ydir','reverse','yticklabel',[])
        set(gca,'xlim',[0 120])
        set(gca,'ytick',0:10:150)
        xlabel('Borehole diameter (mm)')
        title('Fiescherhorn Borehole (Schwerzmann et al 2006)')        
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
end
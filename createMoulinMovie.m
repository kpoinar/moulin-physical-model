%
% for an existing run, create a moulin evolution movie
%
clear variables
%
% Initiate video write
% videofilename = sprintf('Video_Moulin_%s.avi',date);
% vid = VideoWriter(videofilename,'Uncompressed AVI');
% open(vid)
%
%
file2load = 'modeloutputs/02-24-2020/modelrun_R0-2.1_H-500_15d_02-24-2020_1532_outputs.mat';
file2load = 'modeloutputs/02-24-2020/modelrun_R0-2.1_H-500_15d_02-24-2020_1632_outputs.mat';
file2load = 'modeloutputs/02-24-2020/modelrun_R0-2.5_H-500_15d_02-24-2020_1644_outputs.mat';
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-1.9_H-500_12d_02-24-2020_1649_outputs.mat';
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_12d_02-24-2020_1658_outputs.mat';
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_12d_02-24-2020_1708_Froude0.1.mat'; % good movie with ledges emergent
file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_12d_02-24-2020_1738_outputs.mat';
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_12d_02-24-2020_1838_outputs.mat'; % froude 0.001
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_90d_02-24-2020_1841_outputs.mat'; % froude 0.001 and 90 day run
% file2load = 'modeloutputs/02-24-2020/modelrun_R0-0.5_H-500_90d_02-24-2020_1849_outputs.mat'; % froude 1e-9, 90 days


load(file2load);
H = time.z(end);
dlim = 1;  % meters per day x limit for dM, dC, dE plots

% Moulin plot distances
    Moff = -2;  % meters offset to put the moulin more in the center of the plot
    leftlim = -6;
    rightlim = -leftlim;
%
% Colors for moulin plot
watercolor = [0.4 0.8 1];
rockcolor = [0.9 0.7 0.4];
icecolor = [0.93 0.97 1];
pcs = get(0, 'DefaultAxesColorOrder');
%
% Colormaps for time series
Nbrew = 40;
creepcolors = brewermap(Nbrew,'Greens');  % the asterisk means flipud
elastcolors = brewermap(Nbrew,'Blues');
turbucolors = brewermap(Nbrew,'Reds');
%
%
figure(30); clf; set(gcf,'position',[1350  670 1200 675])
% set(gcf,'Color','white')
% four rows
% five columns
%
% 1  2  3  |  4  5
% ---------
% 6  7  8  |  9  10
% 11 12 13 |  14 15
% 16 17 18 |  19 20
%
for ii = length(time.t)%[1:100:length(time.t) length(time.t)]
    %
    
    
    % Melt from turbulence
    h1 = subplot(4,5,[6 11 16],'replace'); hold on
    cc = 0;
    for jj=max(1,ii-Nbrew+2):ii
        cc = cc+1;
        plot((time.dM_major(:,jj)+time.dOC(:,jj))/time.dt*86400,time.z,'-','color',turbucolors(cc,:),'linewidth',2);
    end
    title('Melting by turbulent water')
    set(gca,'xlim',[-dlim dlim]/2)
    xlabel('(meters / day)')
    ylabel('meters above bed')
    set(gca,'xaxislocation','top')
    h1.Position = [h1.Position(1) 0.06 h1.Position(3) 0.54];
    
    
    
    % Opening from creep
    h2 = subplot(4,5,[7 12 17],'replace'); hold on
    cc = 0;
    for jj=max(1,ii-Nbrew+2):ii
        cc = cc+1;
        plot(time.dC_major(:,jj)/time.dt*86400,time.z,'-','color',creepcolors(cc,:),'linewidth',2);
    end
    title('Closing by viscous flow')
    set(gca,'xlim',[-dlim dlim]/2)
    xlabel('(meters / day)')
    %ylabel('meters above bed')
    set(gca,'xaxislocation','top')
    h2.Position = [h2.Position(1) 0.06 h2.Position(3) 0.54];
    
    
    
    
    % Closing from elastic
    h3 = subplot(4,5,[8 13 18],'replace'); hold on
    cc = 0;
    for jj=max(1,ii-Nbrew+2):ii
        cc = cc+1;
        plot(time.dE_major(:,jj)/time.dt*86400,time.z,'-','color',elastcolors(cc,:),'linewidth',2);
    end
    title('Closing by elastic flow')
    set(gca,'xlim',[-dlim dlim]/2)
    xlabel('(meters / day)')
    %ylabel('meters above bed')
    set(gca,'xaxislocation','top')
    h3.Position = [h3.Position(1) 0.06 h3.Position(3) 0.54];
    
    
    
    
    % Qin, Qout, Water level
    h4 = subplot(4,5,1:3,'replace'); hold on
    plot(time.t(1:ii)/86400,time.Qin(1:ii))
%     plot(time.t(1:ii)/86400,time.Qout(1:ii),':')
    plot(time.t(ii)/86400,time.Qin(ii),'.','markersize',24,'color',pcs(1,:))
%     plot(time.t(ii)/86400,time.Qout(ii),'.r','markersize',24)
    set(gca,'xlim',[0 time.t(end)]/86400)
    set(gca,'ylim',[floor(min(time.Qin)) ceil(max(time.Qin))])%[4 8]);
    grid off
    ylabel('Water flow rate (m$^3$/day)')
    xlabel('Days elapsed')    
%   l4 = legend('Water flow into moulin','Water flow out of moulin','location','sw','orientation','horizontal'); l4.Color = 'White';
    l4 = legend('Water flow rate into moulin','location','sw','orientation','horizontal'); l4.Color = 'White';
    h4.Position = [0.075 h4.Position(2) 0.45 h4.Position(4)]; 
    set(h4,'xaxislocation','top')
    %
    h5 = axes; hold on
    set(h5,'position',h4.Position);
    h5.Color = 'None';
    plot(time.t(1:ii)/86400,time.hw(1:ii),'-k')
    plot([0 time.t(ii)],917/1000*H*[1 1],'--k','linewidth',1)
    plot(time.t(ii)/86400,time.hw(ii),'.k','markersize',24)
    set(h5,'yaxislocation','right');
    set(h5,'xlim',[0 time.t(end)]/86400);
    set(h5,'ylim',[0 H])
    ylabel('Water level (m)')
    l5 = legend('Water level','Flotation level','location','se'); l5.Color = 'White';
    grid off
    
    
    % The moulin
    h6 = subplot(4,5,[4 5 9 10 14 15 19 20],'replace'); hold on
    
    % ii=round(length(time.t)/1);
    [~,j] = min(abs(time.z-time.hw(ii)));
    % bed
    patch(leftlim*[-1 -1 1 1 -1],[-100 0 0 -100 -100],rockcolor); hold on
    % water
    patch(Moff+[time.M.xu(1:j,ii); flipud(time.M.xd(1:j,ii)); time.M.xu(1,ii)],[time.z(1:j); flipud(time.z(1:j)); 0],watercolor); hold on
    % ice
    patch(Moff+[time.M.xu(:,ii); leftlim;  leftlim;  time.M.xu(1,ii)], [time.z; time.z(end); time.z(1); time.z(1)],icecolor);
    patch(Moff+[time.M.xd(:,ii); rightlim-Moff; rightlim-Moff; time.M.xd(1,ii)], [time.z; time.z(end); time.z(1); time.z(1)],icecolor);
    %
    % moulin walls
    plot(Moff+time.M.xu(:,ii),time.z,'-k',Moff+time.M.xd(:,ii),time.z,'-k')%,[-1 1]*Mr(j),[1 1]*hw,'-c')
    % ice sheet surface
    plot(Moff+[leftlim time.M.xu(end,ii) NaN time.M.xd(end,ii) rightlim-Moff],[1 1 1 1 1]*H,'-k')
    % ice sheet basal surface
    plot(Moff+[leftlim time.M.xu(1,ii) NaN time.M.xd(1,ii) rightlim-Moff],[0 0 0 0 0],'-k')
    %
    set(gca,'xlim',[leftlim rightlim],'ylim',[-50 H+50]);
    xlabel('Moulin size (m)')
    ylabel('meters above bed')
    set(gca,'yaxislocation','right')
    set(gca,'xaxislocation','top')
    title(sprintf('Day %1.2f',time.t(ii)/86400))
    
    set(gca,'fontsize',30)
    grid off
    h6.Position = [h6.Position(1:3) 0.7];
    
    
    pause(0.01)
    
    frame = getframe(gcf);
    writeVideo(vid,frame);
end


% End video write
close(vid)
% Lauren's plots for the moulin model
function simpleplots(time, save_figures, savelocation,datetime, visible_figures)

%%
timewindow = 1: (43200/time.dt) : (time.sec/time.dt);
%make a few color bars
spect  = brewermap(length(timewindow), 'spectral');
reds   = brewermap(length(timewindow), 'reds');
greens = brewermap(length(timewindow), 'greens');
blues  = brewermap(length(timewindow), 'blues');
oranges  = brewermap(length(timewindow), 'oranges');
greys    = brewermap(length(timewindow), 'Greys');
close 
if visible_figures
    figure
else
    figure('visible', 'off');
end

hold on
set(gcf, 'position', [1          88        2000        1257])
a = subtightplot(1,5,1)

hold on
title('Moulin geometry')
ylabel('Ice thickness')
xlabel('Radius (m)')
axis([-(time.parameters.R0+1) time.parameters.R0+4 0 time.parameters.H])
axis([-4 12 0 time.parameters.H])

for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.M.xd(:,(43200/time.dt).*jj) ,(time.z), 'linewidth', 1, 'color', spect(jj,:) )
plot( time.M.xu(:,(43200/time.dt).*jj) ,(time.z), 'linewidth', 1, 'color', spect(jj,:) )

%plot( time.M.xd(:,(43200/time.dt).*jj) ,(time.z), 'linewidth', 1, 'color', spect(jj,:) )
%plot( time.M.xu(:,(43200/time.dt).*jj) ,(time.z), 'linewidth', 1, 'color', spect(jj,:) )


%plot(zeros(length(timewindow),1) ,    time.hw(:,(43200/time.dt).*jj), 'o', 'markerfacecolor', spect(jj,:), 'markeredgecolor', 'k', 'markersize', 12)
plot( [time.M.xu(round(time.hw(1,(43200/time.dt).*jj)),end), time.M.xd(round(time.hw(1,(43200/time.dt).*jj)),end)], [time.z(round(time.hw(1,(43200/time.dt).*jj))), time.z(round(time.hw(1,(43200/time.dt).*jj)))], '-k','linewidth', 2, 'color',greys(jj,:) )

end
% plot(time.M.r_minor(:,1) ,(time.z), 'linewidth', 3, 'color', 'r')
% plot(time.M.r_minor(:,end) ,(time.z), 'linewidth', 3, 'color', 'r')
% 
% plot(-time.M.r_major(:,1) ,(time.z), 'linewidth', 3, 'color', 'b')
% plot(-time.M.r_major(:,end) ,(time.z), 'linewidth', 3, 'color', 'b')

txt = 'Ice flow \rightarrow';
text(-(time.parameters.R0)+2.5, 0.05*time.parameters.H,txt,'HorizontalAlignment','right', 'fontsize', 20)

plot(time.M.xu(:,1) ,(time.z), 'linewidth', 3, 'color', [0.5 0.5 0.5])


plot(time.M.xd(:,1) ,(time.z), 'linewidth', 3, 'color', [0.5 0.5 0.5])

plot(time.M.xu(:,end) ,(time.z), 'linewidth', 3, 'color', 'k')


plot(time.M.xd(:,end) ,(time.z), 'linewidth', 3, 'color', 'k')




plot( [-(time.parameters.R0+1) time.parameters.R0+12], [time.parameters.H.*0.91, time.parameters.H.*0.91], '--k','linewidth', 1)



b = subtightplot(1,5,2);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (m d^{-1})')
title('Turbulent melting')

for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( ((time.dM_minor(:,(43200/time.dt).*jj) +  time.dM_major(:,(43200/time.dt).*jj))./2)   .* (86400/time.dt) ,(time.z), 'linewidth', 1, 'color', reds(jj,:) )


end





c = subtightplot(1,5,3);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (m d^{-1})')
title('Open channel melting')

for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.dOC(:,(43200/time.dt).*jj)  .* (86400/time.dt)  ,(time.z), 'linewidth', 1, 'color', oranges(jj,:) )


end





d = subtightplot(1,5,4);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (m d^{-1})')
title('Creep closure')
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( ((time.dC_minor(:,(43200/time.dt).*jj) +  time.dC_major(:,(43200/time.dt).*jj))./2)   .* (86400/time.dt)  ,(time.z), 'linewidth', 1, 'color', greens(jj,:) )


end


e = subtightplot(1,5,5);
hold on 
set(gca,'Yticklabel',[])
xlabel('dRadius (m d^{-1})')
title('Elastic deformation')
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( ((time.dE_minor(:,(43200/time.dt).*jj) +  time.dE_major(:,(43200/time.dt).*jj))./2)   .* (86400/time.dt) ,(time.z), 'linewidth', 1, 'color', blues(jj,:) )



set(a, 'position', [0.05 0.1 0.172 0.85])
set(b, 'position', [0.233+0.006 0.1 0.172 0.85])
set(c, 'position', [0.417+0.012 0.1 0.172 0.85])
set(d, 'position', [0.599+0.018 0.1 0.172 0.85])
set(e, 'position', [0.782+0.024 0.1 0.172 0.85])




end



 
if save_figures
  cd(savelocation)
  tmp = datestr(now,'mm-dd-yyyy');
  mkdir(tmp); 
  cd(tmp);
  filename = ['modelrun', '_R0-', num2str(time.parameters.R0), '_H-', num2str(time.parameters.H), '_', num2str(time.parameters.numofdays), 'd_',  datetime, '_geometry.png'];
  saveas(gcf, filename)
end

%%
if visible_figures
    figure
else
    figure('visible', 'off');
end
hold on
set(gcf, 'position', [1          88        1771        1257])
ax(1) = subplot(2,1,1);
hold on
plot(time.t/86400, time.hw) 
plot([time.t(1)/86400 time.t(end)/86400], [time.parameters.H.*0.91 time.parameters.H .*0.91], '--k', 'linewidth',1)
ylabel('Moulin water level (m)', 'fontweight','bold')
xlabel('Time (days)', 'fontweight','bold')

yyaxis right 
plot(time.t/86400, time.S) 
ylabel('Subglacial cross-sectional area (m^2)', 'fontweight','bold')
% 
% 
ax(2) = subplot(2,1,2);
 hold on
 plot(time.t/86400, time.Qin)
 plot(time.t/86400, time.Qout)
 legend('Qin', 'Qout')
 ylabel('Discharge (m^3s^{-1})', 'fontweight','bold')
xlabel('Time (days)', 'fontweight','bold')

linkaxes(ax,'x');

%%

% Plot the moulin

% Moulin plot distances
    Moff = -2;  % meters offset to put the moulin more in the center of the plot
    leftlim = -10;
    rightlim = -leftlim;
%
% Colors for moulin plot
watercolor = [0.4 0.8 1];
rockcolor = [0.9 0.7 0.4];
icecolor = [0.93 0.97 1];

H = time.z(end);

if visible_figures
    figure
else
    figure('visible', 'off');
end

 set(gcf,'position',[377    84   397   721])
hold on

    ii=length(time.t);
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
% 
if save_figures
  cd(savelocation)
  tmp = datestr(now,'mm-dd-yyyy');
  mkdir(tmp); 
  cd(tmp);
  filename = ['modelrun', '_R0-', num2str(time.parameters.R0), '_H-', num2str(time.parameters.H), '_', num2str(time.parameters.numofdays), 'd_',  datetime, '_timeseries.png'];
  saveas(gcf, filename)
end




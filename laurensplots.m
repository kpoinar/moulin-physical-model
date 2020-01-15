% Lauren's plots for the moulin model
function laurensplots(time, save_figures, savelocation,datetime, visible_figures)

%%
timewindow = 1: (43200/time.dt) : (time.sec/time.dt);
%make a few color bars
spect  = brewermap(length(timewindow), 'spectral');
reds   = brewermap(length(timewindow), 'reds');
greens = brewermap(length(timewindow), 'greens');
blues  = brewermap(length(timewindow), 'blues');
close 
if visible_figures
    figure
else
    figure('visible', 'off');
end

hold on
set(gcf, 'position', [1          88        1771        1257])
a = subtightplot(1,4,1)

hold on
title('Moulin radius')
ylabel('Ice thickness')
xlabel('Radius (m)')
axis([0 time.parameters.R0 0 time.parameters.H])
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.Mr(:,(43200/time.dt).*jj) ,(time.z), 'linewidth', 1, 'color', spect(jj,:) )
plot(zeros(length(timewindow),1) ,    time.hw(:,(43200/time.dt).*jj), 'o', 'markerfacecolor', spect(jj,:), 'markeredgecolor', 'k', 'markersize', 12)

end
plot(time.Mr(:,1) ,(time.z), 'linewidth', 3, 'color', 'k')

plot(time.Mr(:,end) ,(time.z), 'linewidth', 3, 'color', 'k')
plot( [0, time.parameters.R0], [time.parameters.H.*0.91, time.parameters.H.*0.91], '--k','linewidth', 1)
b = subtightplot(1,4,2);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (m d^{-1})')
title('Turbulent Melting')

for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.dM(:,(43200/time.dt).*jj) .* (86400/time.dt) ,(time.z), 'linewidth', 1, 'color', reds(jj,:) )


end



c = subtightplot(1,4,3);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (m d^{-1})')
title('Creep closure')
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.dC(:,(43200/time.dt).*jj) .* (86400/time.dt),(time.z), 'linewidth', 1, 'color', greens(jj,:) )


end


d = subtightplot(1,4,4);
hold on 
set(gca,'Yticklabel',[])
xlabel('dRadius (m d^{-1})')
title('Elastic deformation')
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.dE(:,(43200/time.dt).*jj) .* (86400/time.dt) ,(time.z), 'linewidth', 1, 'color', blues(jj,:) )


end



set(a, 'position', [0.08 0.1 0.2 0.8])
set(b, 'position', [0.3 0.1 0.2 0.8])
set(c, 'position', [0.52 0.1 0.2 0.8])
set(d, 'position', [0.74 0.1 0.2 0.8])
 
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
subplot(2,1,1)
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
 subplot(2,1,2)
 hold on
 plot(time.t/86400, time.Qin)
 plot(time.t/86400, time.Qout)
 legend('Qin', 'Qout')
 ylabel('Discharge (m^3s^{-1})', 'fontweight','bold')
xlabel('Time (days)', 'fontweight','bold')


% 
if save_figures
  cd(savelocation)
  tmp = datestr(now,'mm-dd-yyyy');
  mkdir(tmp); 
  cd(tmp);
  filename = ['modelrun', '_R0-', num2str(time.parameters.R0), '_H-', num2str(time.parameters.H), '_', num2str(time.parameters.numofdays), 'd_',  datetime, '_timeseries.png'];
  saveas(gcf, filename)
end




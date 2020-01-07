% Lauren's plots for the moulin model
function laurensplots(time)
%%
timewindow = 1: (43200/time.dt) : (time.sec/time.dt);
%make a few color bars
spect  = brewermap(length(timewindow), 'spectral');
reds   = brewermap(length(timewindow), 'reds');
greens = brewermap(length(timewindow), 'greens');
blues  = brewermap(length(timewindow), 'blues');

figure
hold on
set(gcf, 'position', [1          88        1771        1257])
a = subtightplot(1,4,1)

hold on
title('Moulin radius')
ylabel('Ice thickness')
xlabel('Radius (m)')
axis([0 time.parameters.R0 0 time.parameters.H])
for jj = 1:1:length(timewindow) % plot every 12 hours for the full model run
plot( time.Mr(:,(43200/time.dt).*jj),time.z, 'linewidth', 1, 'color', spect(jj,:) )


end

b = subtightplot(1,4,2);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (mm h^{-1})')



c = subtightplot(1,4,3);
hold on 
set(gca,'Yticklabel',[]) 
xlabel('dRadius (mm h^{-1})')

d = subtightplot(1,4,4);
hold on 
set(gca,'Yticklabel',[])
xlabel('dRadius (mm h^{-1})')



set(a, 'position', [0.08 0.1 0.2 0.8])
set(b, 'position', [0.3 0.1 0.2 0.8])
set(c, 'position', [0.52 0.1 0.2 0.8])
set(d, 'position', [0.74 0.1 0.2 0.8])
 



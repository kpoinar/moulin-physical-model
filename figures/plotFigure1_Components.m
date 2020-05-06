% Create Figure 1, which illustrates the multi processes that determine the
% geometry of a moulin
%
% kristin poinar
% march 21, 2020
% shelter in place
%
clear variables
%
% Load moulin file
load modeloutputs/03-21-2020/modelrun_R0-1_H-800_1d_03-21-2020_1752_outputs.mat
% A short (1 day) run, has time struct with fields
% dC_major   CREEP
% dC_minor   CREEP
% dM         TURB MELTING
% dOC        WATERFALL MELTING
% dE_major   ELASTIC
% dE_minor   ELASTIC
% dG         GLENS
% dP         NOT USED
% It is missing dF, from freezing (not calculated)
%%
% Some colors
watercolor = [0.1 0.8 1];
blueicecolor = [0.1 0.5 1];
whiteicecolor = [0.92 0.95 1];
rock = [0.6824 0.4667 0];
%
% Define the time step to plot.  It took the subglacial system 0.05 days 
% to recover from initial conditions.
ii = find(time.t > 0.05*86400,1,'first');
% Move it forward a little
ii = ii+37;
%
wet = locatewater(time.hw(ii),time.z);
H = time.z(end);
    
figure(10); clf; hold on
set(gcf,'position',[194  169  1151  636]);

wall0 = -time.M.r_minor(:,ii);
wall1 = -wall0;
%
% Which processes are going to be shown?
changes = [time.dC_minor(:,ii)  time.dE_minor(:,ii)  time.dM(:,ii)  time.dOC(:,ii)  time.dG(:,ii)];
% Maximum change is 0.003 meters.  Scale this up to an avg change of 0.5 m:
avgchg = 0.25;
changescaled = changes .* repmat(avgchg./mean(changes(wet,:)),numel(time.z),1) .* sign(mean(changes));
% fix the open channel part
changescaled(:,4) = avgchg./mean(changes(~wet,4)) .* changes(:,4);
%
% process labels to title the moulins
labs = {'Viscous','Elastic','Turbulence','Open Channel','Glens Law','New Geometry'};
%
% Create x offsets to make 5 moulins that are offset from one another
doff = 4;
offset = 0 : doff : 25;
% Plot bedrock
patch([-2 22 22 -2 -2],[-50 -50 0 0 -50],rock)
grid off
start = -2;
% plot colors
pcs = [0.0706    0.1961    0.9   ; ...
       0.8       0.2       0.2001; ...
       0.7451    0.7451    0.1;    ...
       0.8510    0.3725    0.0078; ...
            0    0.4902    0.1961; ...
       0.9       0.2       0.4];
%
% Lay the background ice color
patch([-2; wall0; -2; -2],[0; time.z; H; 0],whiteicecolor)
for pp = 1:4%5    
    offs = offset(pp+1);
    patch([wall0 + offs; flipud(wall1) + offs-doff; wall0(1)],[time.z; flipud(time.z); 0],whiteicecolor)
end
         


% Finally, add a "total change" moulin that shows the contribution from all
% sources...
pp = 6;

offs = offset(pp);
% Ice surface
%plot([start offs+wall0(end)],[H H],'-k')
start = wall1(end)+offs;
%
% Calculate the change to the left moulin wall from all processes
scale = 250;
gscl = 1.0;  % reduce GFL 
chgL = (-sum(changes(:,1:4),2) + gscl*changes(:,5)) * scale;   
% Calculate the change to the right moulin wall from all processes
% (remembering NOT to add open channel change, #4, to this wall)
chgR = (sum(changes(:,1:3),2) + gscl*changes(:,5)) * scale;
% Water (with walls after all processes)
patch(offs+[wall0(wet)+chgL(wet); flipud(wall1(wet)+chgR(wet)); wall0(wet(1))],[time.z(wet); flipud(time.z(wet)); time.z(wet(1))],watercolor);

% Moulin walls (after all processes) 
% Left
plot(offs+wall0+chgL+0*chgL(end),time.z,'-k','linewidth',3)%,'color',pcs(pp,:))
% Right
plot(offs+wall1+chgR-0*chgR(end),time.z,'-k','linewidth',3)%,'color','pcs(pp,:))
%
% Ice surface
%plot([offs+wall1(end) offs+2],[H H],'-k')


% White ice background
% right of moulin 5
patch(0*doff+[wall1+offs-doff; flipud(wall0+chgL)+offs; wall1(1)+offs-doff],[time.z; flipud(time.z); 0],whiteicecolor)
% right of moulin 6
patch([22; 22; flipud(wall1+chgR)+offs; 22],[0; H; flipud(time.z); 0],whiteicecolor)

% Moulin walls (before all processes)
plot(offs+[wall0 wall1],time.z,'--k','linewidth',2)
%




% Now return to the 1st panel and plot the "starter" moulin (black), the
% "starter" water (blue), and the changes to it (colorful new walls)

% Step through all processes and plot the moulin, water, and the change
% to both moulin walls from each process
for pp = 1:5
    offs = offset(pp);
    % Water 
    patch(offs+[wall0(wet); flipud(wall1(wet)); wall0(wet(1))],[time.z(wet); flipud(time.z(wet)); time.z(wet(1))],watercolor);
    % Moulin walls (beginning of timestep)
    plot(offs+[wall0 wall1],time.z,'-k')
    % Ice surface
    %plot([start offs+wall0(end)],[H H],'-k')
    start = wall1(end)+offs;
    %
    % Now plot the change to the moulin from this process
    chg = changescaled(:,pp);
    %chg(chg==0) = NaN;
    %
    % Left wall
    if pp<5 % then it's radius change, I subtract to get left wall position
        plot(offs+wall0-chg,time.z,'-','color',pcs(pp,:))
    else % then it's GFL, I add to get left wall position
        % and make it a little larger
        chg = chg*2;
        plot(offs+wall0+chg,time.z,'-','color',pcs(pp,:))
    end
    % Right wall
    if pp == 4, 
        % Open channel only affects the upstream wall
        plot(offs+wall1+0*chg,time.z,'-','color',pcs(pp,:))
    else
        plot(offs+wall1+chg,time.z,'-','color',pcs(pp,:))
    end
    
    text(offs,H+30,labs{pp},'color',pcs(pp,:),'fontweight','bold','fontsize',20,'horizontalalignment','center')
end
plot([0 1]+start,[H H],'-k')


set(gca,'ylim',[-50 H+50])
set(gca,'xlim',[-2 22],'ylim',[-50 H+70])
set(gca,'xtick',[0 : 2 : 24])
xlabel('Horizontal distance (m)')
ylabel('Height above bed (m)')

% Total change label
pp = 6;
text(offset(pp),H+30,labs{pp},'color','k','fontweight','bold','fontsize',20,'horizontalalignment','center')

box off
%%
figfile = sprintf('Figure1_ModelComponents_%s',date);
savefig(strcat(figfile,'.fig'))
print(figfile,'-dtiff')



% Some notes after the first version of this script:
% On the right-most panel, it'd be cool to illustrate the water as it would 
% be in the moulin "after" all the geometry changes.
% In the current case, the water would be higher.
% However, this isn't what our model actually does.  It sets the water
% level according to the Schoof model outputs.
% How to illustrate???
% Could it be as simple as drawing the new walls and water at the current
% level?
            
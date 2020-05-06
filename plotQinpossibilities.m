% Plot Qin range, Qin base possibilities versus elevation
% for comparison with Lauren's reanalysis melt data by elevation
%
function plotQinpossibilities(H,Qinrange,Qinbase,success)%,fig2open)

fig2open = 'possibleQin.fig';
figpath = './modeloutputs/';

openfig(strcat(figpath,fig2open));

hold on

for ii=1:length(success)
    
    switch success(ii)
        case 1, color = 'b';
        case 0, color = 'r';
    end
    
    n=3;
   plot(H*[1 1]-length(success)/2*n + n*ii,max([0 0],Qinrange(ii)*[-0.5 0.5] + Qinbase(ii)),color) 
    
end


xlabel('Ice thickness (m)')
ylabel('Qin (m$^3$/s)')

savefig(strcat(figpath,fig2open));



if 0
    
   Qmeantable = [3.7 3.2 4.5 4.5 4.2 4.6 3.2 3.4 3.4 3.7];
   Qrangetable = [5 2.3 3.2 2.1 1.7 1.1 0.8 0.7 0.7 0.8];
   Htable = [669 820 947 1058 1159 1252 1339 1420 1497 1569];
   
   for ii=1:length(Htable)
       plot(Htable(ii)*[1 1],Qmeantable(ii)+[-0.5 0.5]*Qrangetable(ii),'.-k','markersize',40)
   end
    
end
%
% Hansen & Landauer 1958
%
clear all
%{
% Depth (feet and meters)
Zfeet = [185:50:435 535:50:985 1010 1035:50:1135];
Z = Zfeet / 2 * 12 * 2.54 / 100;
%
% Pressure (kg/cm^2 and Pa)
Pcgs = [3.58 5.30 6.70 8.10 9.50 10.90 13.70 15.10 16.50 17.90 19.31 ...
    20.71 22.11 23.55 24.92 26.32 27.02 27.72 29.13 30.53];
%
% Hole diameter (inches and cm)
Dinch = [5.836 5.824 5.780 5.808 5.776 5.758 5.934 5.518 5.470 5.372 ...
    5.140 5.000 4.958 4.474 4.410 4.022 3.786 3.684 3.068 2.704];
R = Dinch / 2.54 / 2;
%
% Strain
eps = (R-R(1)) / R(1);

figure(100); 
subplot(1,2,1,'replace'); plot(R,Z); set(gca,'ydir','reverse')
subplot(1,2,2,'replace'); plot(eps,Z); set(gca,'ydir','reverse')

%
%
%}
% Meighen Ice Cap
G = makeGlobalParams;
t = [5 4 2 1]*G.seconds;
Z = [108; 93; 78; 48];
epsmat = [29 19.5 7.3 3; 13.5 9.7 3.1 2; 7 5 2 1; 2.5 2 1.5 1];
a0 = 1;
amat = a0 * exp(-epsmat);

figure(101); 
subplot(1,3,1,'replace'); plot(t/G.seconds,amat,'--o')
    xlabel('Time (years)')
    ylabel('Borehole radius / original radius')
    legend('108 m','93 m','78 m','48 m')
subplot(1,3,2,'replace'); plot(amat,Z,'-o'); set(gca,'ydir','reverse')
    xlabel('Borehole radius / original radius')
    ylabel('Depth (m)')
    legend('5 years','4 years','2 years','1 year','location','se')
    
    
sigZ = -G.rhoi*G.g*Z;
subplot(1,3,3,'replace'); plot(epsmat,sigZ/1e3,'o-')
    xlabel('Strain')
    ylabel('\sigma_H (kPa)')
    legend('5 years','4 years','2 years','1 year','location','ne')
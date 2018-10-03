function C = makeConstants()

C.To = 273.15; % Kelvin
C.rhoi = 910;  % kg/m3
C.rhow = 1000; % kg/m3
C.ki = 2.1;    % J/mKs
C.cp = 2115;   % J/kgK
C.Lf = 335000; % J/kg
C.g = 9.8;     % m/s2
C.E = 1e9;     % Pa; Young's elastic modulus (Vaughan 1995)






%C.E = 1e12;







C.nu = 0.3;    % []; Poissons ratio for ice
%
% Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35
% A = Astar * exp(-Q / R *(1/Th - 1/Tstar))
% Th = 263 + 7e-8*P
% Tstar = T + 7e-8*P
C.Tstar = 263;
C.Astar = 3.5e-25;
C.c = 7e-8;
C.Qcless = 6e4;
C.Qcmore = 11.5e4;    
C.R = 8.314;                % ideal gas constant
C.a = 7e-8;
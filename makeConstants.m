function C = makeConstants()

C.To = 273.15; % Kelvin
C.rhoi = 910;  % kg/m3; Ice density
C.rhow = 1000; % kg/m3; Water density
C.ki = 2.1;    % J/mKs
C.cp = 2115;   % J/kgK
C.Lf = 335000; % J/kg; Latent heat of fusion
C.g = 9.8;     % m/s2; Gravity
C.E = 1e9;     % Pa; Young's elastic modulus (Vaughan 1995)
C.A = 6e-24;   % 1/Pa3/s; Glen's law fluidity coefficient (Schoof 2010)
C.f = 0.1;     % unitless; Darcy-Weisbach friction factor (0.1 in Matt's code, 0.0375 in Schoof 2010)
C.n = 3        % unitless; Glen's law exponent (Schoof 2010)


% Subglacialsc model constants
C.c1 = 1/pi/C.Lf % units; Melt opening parameter (Schoof 2010)
C.c2 = C.A*C.n^(-C.n) % units; Closure parameter (Schoof 2010)
C.c3 = 2^(1/4) * (pi+2)^(1/2)/(pi^(1/4) * (C.rhow*C.f)^(1/2)) % units; Flux parameter (Schoof 2010)







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
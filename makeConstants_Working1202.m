function C = makeConstants()
 
C.T0 = 273.15; % Kelvin
C.rhoi = 910;  % kg/m3; Ice density
C.rhow = 1000; % kg/m3; Water density
C.ki = 2.1;    % J/mKs
C.cp = 2115;   % J/kgK
C.Lf = 335000; % J/kg; Latent heat of fusion
C.g = 9.8;     % m/s2; Gravity
C.E = 1e9;     % Pa; Young's elastic modulus (Vaughan 1995)
C.A = (8e-24);   % 1/Pa3/s; 6e-24 Glen's law fluidity coefficient (Schoof 2010)
C.f = 0.1;     % unitless; Darcy-Weisbach friction factor (0.1 in Matt's code, 0.0375 in Schoof 2010)
C.n = 3;        % unitless; Glen's law exponent (Schoof 2010)


% Subglacialsc model constants
C.c1 = 1/C.rhoi/C.Lf; % units; Melt opening parameter (Schoof 2010)
C.c2 = 2*C.A*C.n^(-C.n); % units; Closure parameter (Schoof 2010)
%C.c3 = 2^(1/4) * (pi+2)^(1/2) / (pi^(1/4) * (C.rhow*C.f)^(1/2)); % units; Flux parameter (Schoof 2010)
%      2^(5/4) / pi**(1/4) * sqrt( pi/ ((pi+2)*rho_w*f) ) # from Matt Covington pdf, 
% for a semi-circular conduit, modified from Schoof whose equation appears to be incorrect. 
C.c3 = ( 2.^(5./4) ./pi.^(1./4.) .* sqrt(pi./(pi + 2)) )/ sqrt(C.rhow*C.f);


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

% Turbulence parameters
C.mu = 0.0017916; % Pa*s at 0∞C, https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
C.kw = 0.555; % J/(mKs) at 0.01∞C, https://www.engineeringtoolbox.com/water-liquid-gas-thermal-conductivity-temperature-pressure-d_2012.html
C.cw = 4210  ; %J / (kg * K)   heat capacity of water for unit mass,  from Jarosch & Gundmundsson (2012)

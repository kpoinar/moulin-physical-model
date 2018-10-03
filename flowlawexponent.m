function A = flowlawexponent(T,P,C)
%
% Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35
% A = Astar * exp(-Q / R *(1/Th - 1/Tstar))
% Th = 263 + 7e-8*P
% Tstar = T + 7e-8*P
%
% Takes temperature, pressure, and a struct containing these constants (C):
Tfrac = 1./(T + C.a*P) - 1./(C.Tstar + C.a*P);
Qc = C.Qcless * ones(size(T));
Qc(T>C.Tstar) = C.Qcmore;
A = C.Astar * exp(-Qc/C.R .* Tfrac);
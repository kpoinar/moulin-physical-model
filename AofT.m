function A = AofT(T)
% Calculate the flow law parameter, A, given ice temperature T
% A is in Pa^-3 s^-1
% T must be in degrees Celsius
%
Alookup = [2.4 1.7 0.93 0.35 0.21 0.12 0.037 0.01 0.0026]*1e-24; % s-1 Pa-3 (Cuffey Table B.1)
Tlookup = -[0 2 5 10 15 20 30 40 50]; % °C

A = interp1(Tlookup,Alookup,T,'pchip');
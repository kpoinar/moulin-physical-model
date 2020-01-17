% Creep closure of a water-filled borehole
%
% Based on boreholeclosure/HomeworkProblem_Vostok3G.m
% which I did in 2013 for crevasse model
%
% Borehole 3G at Vostok
% by Blinov and Dmitriev (1987) and Salamatin et al (1998)
% From Table 4 in Talalay and Hooke, 2007 (Annals)
%
function dC = creep(Mr,z,H,stress,T,dt,E,C)
%
% Assume A for temperate ice for now
%A = mean(G.A0less * exp(-G.Qless / G.R ./ (G.To-[0 15])));
%A = G.A0less * exp(-G.Qless / G.R ./ (T));
% Pressure:
P = C.rhoi*C.g*(H-z);
A = flowlawexponent(mean(T,2),P,C) ;

% % Ice hydrostatic stress (closure)
% sigzi = C.rhoi*C.g*(H-z);
% % Water hydrostatic stress (opening)
% sigzw = -C.rhow*C.g*(hw-z);
% sigzw(~wet) = 0; % Anywhere that is not wet does not have the opening force from water
% ^^ The above is now calculated in the main function, moulingeom.m, and
% passed to this function
% In the struct sig,
% water pressure is positive (opening)
% ice pressure is negative (closing)
%
% Total stress
sigmaZ = stress.cryo + stress.hydro;
epsdot = E*A .* (sigmaZ/3).^3;  % boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
%
% figure(6); 
% subplot(1,2,1); plot(sigmaZ,z);
% subplot(1,2,2); plot(epsdot,z);
%
% Borehole closure
Mrnew = Mr .* exp(epsdot * dt);
%Mrnew = max(Mrnew,0);
% Creep closure rate
dC = Mrnew-Mr;
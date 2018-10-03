% Creep closure of a water-filled borehole
%
% Based on boreholeclosure/HomeworkProblem_Vostok3G.m
% which I did in 2013 for crevasse model
%
% Borehole 3G at Vostok
% by Blinov and Dmitriev (1987) and Salamatin et al (1998)
% From Table 4 in Talalay and Hooke, 2007 (Annals)
%
function dC = creep(Mr,z,H,hw,T,dt,E,C)
%
% Assume A for temperate ice for now
%A = mean(G.A0less * exp(-G.Qless / G.R ./ (G.To-[0 15])));
%A = G.A0less * exp(-G.Qless / G.R ./ (T));
% Pressure:
P = C.rhoi*C.g*(H-z);
A = flowlawexponent(mean(T,2),P,C);

% Ice hydrostatic stress (closure)
sigzi = C.rhoi*C.g*(H-z);
% Water hydrostatic stress (opening)
sigzw = -C.rhow*C.g*(hw-z);
sigzw(z>hw) = 0;
% Total
sigmaZ = sigzi + sigzw;
epsdot = E*A/1 .* (sigmaZ/3).^3;  % boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
%
% figure(6); 
% subplot(1,2,1); plot(sigmaZ,z);
% subplot(1,2,2); plot(epsdot,z);
%
% Borehole closure
Mrnew = Mr .* exp(-epsdot * dt);
%Mrnew = max(Mrnew,0);
% Creep closure rate
dC = Mrnew-Mr;
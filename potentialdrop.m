function dP = potentialdrop(Qin,wet,eggperi,dt,C,f)
%
% Convert potential energy released into latent heat to melt the ice that
% is NOT underwater
% dP = C.rhow / C.rhoi *C.g / C.Lf * Qin*dt ./ (2*pi*Mr);
%
% This is a very large amount of melting.  In fact, NOT all of the PE will
% transfer to the walls; the water will quickly reach terminal velocity (so
% only some is lost to KE which transfers to the water surface on impact).
% From there, the water will frictionally heat the air, and some fraction
% of that heat will then reach the moulin walls.  The fraction is unknown,
% but should scale inversely with moulin radius: the farther away the wall,
% the less likely the air will fully transfer the energy.
% For now, try a fraction of energy transfer f:
% f ~ 1/r
%
%f = 1;%./Mr;
%Mr = ones(length(Mr));
dP = C.rhow / C.rhoi *C.g / C.Lf * Qin*dt ./ (eggperi);

%dP = C.rhow / C.rhoi *C.g / C.Lf * Qin*dt ./ (2*pi*);


dP = dP .* f;
%
% Erase dP below the water line
dP(wet) = 0;
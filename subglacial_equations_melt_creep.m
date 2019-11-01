function f = subglacial_equations_melt_creep(t,R,hw,S)
%IcePressure

%Moulin cross-section area in function of z
A_R = pi*r^2;
%Head partial differential equation
dhdt = 1/A_R * ( R - C3*S(5/4) * sqrt((rho_w*g*h)/L) );
%Channel cross section area partial differential equation
dSdt = C1 * C3 * S^(5/4) * ((rho_w*g*h)/L)^(3/2) - C2 * ( Pi - rho_w*g*h )^n * S;
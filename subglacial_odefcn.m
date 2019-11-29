function dydt = subglacial_odefcn(t,y,Mr,z,Qin,H,dx,C)

%{ 
This function is based on Schoof 2010, without the cavity. 
This function is the input of the ODE solver for the subglacial system .

Inputs:
-------
t = time (s)
y = needs to be there for the solver
Mr = [vector] moulin radius in function of z (m)
z = vertical position of modeled Moulin radius (m)
Qin = Recharge in funtion of time or not (m3/s)
H = [single variable] Ice thickness (m)
dx = [single variable] distance between terminus and moulin. Equivalent to channel length (L in schoof)? (m)

Outputs:
--------
dydt(1): dhdt. Moulin head in function of time
dydt(2): dSdt. Channel cross-section area in function of time

%}
hw = y(1);
S = y(2);
%y(1) = hw; I don't know if doing it this way would work
%y(2) = S; 

%Import constants from the matlab file makeConstants.
%C = makeConstants;
%Create empty dydt (necessary for it to work)
dydt = zeros(2,1); 
%WaterPressure
Pw = C.rhow .* C.g .* hw;
%IcePressure 
Pi = C.rhoi .* C.g .* H;
%Moulin cross-section area in function of z
Mr_z = interp1(z,Mr,hw);
A_R = pi .* Mr_z.^2;
%disp(Mr_z)
%Recharge in moulin in function of time
%
%Head partial differential equation --- dhdt
dydt(1) = 1./A_R .* ( Qin - C.c3.*S.^(5/4) .* sqrt(Pw./dx) );
%Channel cross section area partial differential equation --- dSdt
dydt(2) = C.c1 .* C.c3 * S.^(5./4) .* (Pw./dx).^(3./2) ...
         - C.c2 .* ( Pi - Pw ).^C.n .* S;

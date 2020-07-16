%function [t,hw,S,Qout] = subglacialsc(Mr,z,Qin,H,dx,C,tspan,y0) %for
%TestSubglacial.m
function [hw,S,Qout, dydt_out] = subglacialsc(Ms,z,Qin,Qvadd,H,L,C,dt, tspan,y0, opt)
% import constants
%C = makeConstants;

%{
This function solves the coupled differential equations for the subglacial part
It needs the function subglacial_equations_melt_creep to work.


input:
------
Mr = Moulin radius in function of z. 
z = Vertical position of modeled Moulin radius (m)
Qin = Supraglacial recharge. Can be constant or a function of time.
H = Ice thickness (m)
dx = Distance between terminus and moulin. Equivalent to channel length L in schoof? (m)
tspan = [t0, tf] vector composed of t initial and t final. 
Which will be the current timestep in the single ode solver, and the next timestep
y0 = initial condition
    y0(1) = initial/previous timestep value of hw (moulin head) 
    y0(2) = initial/previous timestep value of S (cross-section area) or the 

output:
-------
hw = Moulin head (m)
S = Channel cross-section area (semi-circular)
Qout = Discharge out of the channel


We use ode23s to solve the ODE. is it a stiff ode? 
With Matt Covington, we've been using LSODA in python, and ode23s seems to
be the equivalent in matlab. 
-look into that

%}

%[t,y] = ode45(@(t,y) subglacial_odefcn(t,y,Mr,z,Qin,H,dx,C), tspan, y0);

Qin = Qin+Qvadd; %adds or remove the volume of water squeezed or relaxed when the moulin creep or elastic.

%Using ode15s deals with the equation stiffness problem

[t,y] = ode15s(@(t,y) subglacial_odefcn(t,y,Ms,z,Qin,H,L,C), tspan, y0, opt); %testing to see if a stiff solver will deal with the timesteping issues

%tic, ode15s(@(t,y) subglacial_odefcn(t,y,Ms,z,Qin,H,L,C), tspan, y0, opt), toc %testing to see if a stiff solver will deal with the timesteping issues


hw = y(end,1); % moulin head (m)
S = y(end,2); % channel cross-section area
dydt_out = NaN;

if hw == y0(1)



Pw =  C.rhow .* C.g .* hw;
Pi = C.rhoi .* C.g .* H;
    
dydt_out = C.c1 .* C.c3 * S.^(5./4) .* (Pw./L).^(3./2) ...
         - C.c2 .* ( Pi - Pw ).^C.n .* S;


 S = S + dydt_out*dt;

end

Qout = C.c3 .* S.^(5/4) .* sqrt( C.rhow*C.g*hw/L); % discharge out the channel



% Elastic deformation of a cylindrical hole in a material, with water
% pressure Pw and ice pressure Pi
% Based on Aadnoy 1987: Model for Fluid-Induced and In-Situ Generated
% Stresses in a Borehole (in rock)
% 
% Created July 19, 2018 by Kristin Poinar for use in the moulin model 
% 
% Modified November 1, 2019 to fix the many errors in the equation I
% derived from Aadnoy.
%
% Updated April 2, 2022 to fix the error that it is the CHANGE in pressure,
% not the raw pressure, that drives the radial deformation dr
%
%
% This solution assumes plane strain at the base of the ice sheet (0 vertical strain; Aadnoy assumes
% the necessary vertical stress to make that true)
%
% Inputs:   Mr      Moulin radius (m)
%           stress  Struct containing stress components at this timestep (Pa) 
%           stress_prev                     ... at the previous timestep (Pa)
%           C       Struct containing physical constants
% 
% Outputs:  dr      Change in moulin radius over this timestep, due to the
%                   change in stress (m)


function dr = elastic(Mr,stress,stress_prev,C)

    % Water pressure ("outward")
    %Pw = C.rhow * C.g * (hw-z); % Pa
    % Ice pressure ("inward")
    %Pi = C.rhoi * C.g * (H-z);  % Pa
    % Total pressure (inward, unless water level is above flotation)
    %P = Pw - Pi;  % Pa
    
    
    dP = stress.hydro + stress.cryo - (stress_prev.hydro + stress_prev.cryo);
    %
    % Radial deformation (relax to equilibrium)
    % Without a/E prefactor:
    dr = (1 + C.nu)*(dP - 0.5*(stress.sigx+stress.sigy)) + 0.25 * (stress.sigx-stress.sigy)*(1 - 3*C.nu - 4*C.nu^2) + 0.25 * stress.tauxy * (2 - 3*C.nu - 8*C.nu^2);
    % Radial deformation (moulin radius a or Mr):
    dr = dr .* Mr / C.Y;  % meters
    

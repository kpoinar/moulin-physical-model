% Elastic deformation of a cylindrical hole in a material, with water
% pressure Pw and ice pressure Pi
% Based on Aadnoy 1987: Model for Fluid-Induced and In-Situ Generated
% Stresses in a Borehole (in rock)
% 
% Created July 19, 2018 by Kristin Poinar for use in the moulin model 
% Modified November 1, 2019 to fix the many errors in the equation I
% derived from Aadnoy.
%
%
% This solution assumes plane strain at the base of the ice sheet (0 vertical strain; Aadnoy assumes
% the necessary vertical stress to make that true)
%
function dr = elastic(Mr,stress,C)

    % Water pressure ("outward")
    %Pw = C.rhow * C.g * (hw-z); % Pa
    % Ice pressure ("inward")
    %Pi = C.rhoi * C.g * (H-z);  % Pa
    % Total pressure (inward, unless water level is above flotation)
    %P = Pw - Pi;  % Pa
    P = stress.hydro + stress.cryo;
    %
    % Radial deformation (relax to equilibrium)
    % Without a/E prefactor:
    dr = (1 + C.nu)*(P - 0.5*(stress.sigx+stress.sigy)) + 0.25 * (stress.sigx-stress.sigy)*(1 - 3*C.nu - 4*C.nu^2) + 0.25 * stress.tauxy * (2 - 3*C.nu - 8*C.nu^2);
    % Radial deformation (moulin radius a or Mr):
    dr = dr .* Mr / C.E;  % meters
    

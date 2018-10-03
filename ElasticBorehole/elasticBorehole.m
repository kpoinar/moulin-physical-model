% Elastic deformation of a cylindrical hole in a material, with water pressure Pw
% Based on Aadnoy 1987: Model for Fluid-Induced and In-Situ Generated
% Stresses in a Borehole (in rock)
%
%
% plane strain at the base of the ice sheet (0 vertical strain; necessary
% vertical stress to make that true)
%
clear all
G = makeGlobalParams;
sigx = 100e3;
sigy = -100e3;
nu = 0.3;
tauxy = 100e3;
Pw = G.rhow*G.g*500;  % 500 m of water pressure
a = 1;  % borehole radius
E = 1e9; % Young's elastic modulus (Vaughan 1995)
%
%ni = 100;
%avec = zeros(1,ni);
%for i=1:ni
    % Radial deformation (relax to equilibrium)
    dr = (sigx-sigy)*((3-nu)/4*a-a^2*nu) + (sigx+sigy)*a/2*(1+nu) + Pw*(nu-0/5)*a + tauxy*a*(3/4 - nu/2 - 2*nu^2);
    dr = dr / E;
    %avec(i) = a;
    %a = a - dr;
%end
fprintf('dr = %1.1e m of deformation in borehole of radius a = %1.1e m\n',dr,a)
%figure(1); clf; plot(avec)
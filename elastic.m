% Elastic deformation of a cylindrical hole in a material, with water pressure Pw
% Based on Aadnoy 1987: Model for Fluid-Induced and In-Situ Generated
% Stresses in a Borehole (in rock)
%
%
% plane strain at the base of the ice sheet (0 vertical strain; necessary
% vertical stress to make that true)
%
function dr = elastic(sigx,sigy,nu,E,tauxy,Pw,Mr)

    % Radial deformation (relax to equilibrium)
    dr = (sigx-sigy)*((3-nu)/4*Mr-Mr.^2*nu) + (sigx+sigy)*Mr/2*(1+nu) + Pw*(nu-1)*Mr + tauxy*Mr*(3/4 - nu/2 - 2*nu^2);
    dr = dr / E;
    %
    % The deformation acts to close the borehole:
    %dr = -dr;
    % Yes, in fact, it does; so dr is negative.
    % I had had a typo in the above long equation:
    % correct: ... + Pw*(nu-0.5)*Mr + ...
    % typo   : ... + Pw*(nu-0/5)*Mr + ...
    % which made this large term flip sign (correct is +, was showing -)
    % so I adjusted the sign of dr when it came out dr > 0
    %
    Mrnew = max(0,Mr + dr);
    dr = Mrnew - Mr;
    %
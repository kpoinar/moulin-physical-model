function dV = calculate_dQ_deformation(dMr_major,dMr_minor,M,z,wet)
	% Calculate the change of volume per second produce by moulin deformation.
	%dMr_major: dE_major or dC_major
	%dMr_minor: dC_minor or dC_minor

	%Old ellipse area - New ellipse area
    dA_ellipse = pi .* M.r_major .* M.r_minor - pi .* (M.r_major+dMr_major) .* (M.r_minor+dMr_minor);
    %Old circle area - New circle area
    dA_circle = pi .* M.r_minor.^2 - pi .* (M.r_minor+dMr_minor).^2;
    %total change in area calculated with half the ellipse and half the circle
    dA_tot = 0.5.*dA_ellipse + 0.5.*dA_circle;
    %change in volume below the water level for an entire timestep!
    dV = trapz(z(wet),dA_tot(wet));
 
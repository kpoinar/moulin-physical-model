%create a new elastic module based on revised equation 7



function dr = elastic2(Mr, stress, hw_einit, Mr_einit, dt ,C, z, wet)

tmpstresshydro = C.g .* C.rhow .* (hw_einit - z) ;
tmpstresshydro(~wet) = 0;

P0 = tmpstresshydro + stress.cryo;

P1 = stress.hydro + stress.cryo;

dP_dt =  ((P1 - P0) / dt);

dMr_dt = ((Mr - Mr_einit) / dt); 


%dr = ((1 + C.nu) ./ C.Y ) .* ((Mr .* dP_dt) + (P1 .* dMr_dt)) .* dt ;

dr = (Mr .* (1 + C.nu) .* dP_dt + ((1 + C.nu) .* (P1 - 0.5*(stress.sigx+stress.sigy)) + 0.25 * (stress.sigx-stress.sigy)*(1 - 3*C.nu - 4*C.nu^2) ...
     + 0.25 * stress.tauxy * (2 - 3*C.nu - 8*C.nu^2)) .* dMr_dt);

dr = dr * (dt/C.Y);



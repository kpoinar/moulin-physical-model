% dE_minor = res{1,2}.dE_minor;
% dE_major = res{1,2}.dE_major;
% z = res{1,2}.z;
% time = res{1,2}.t;
% 
% 
% figure(202)
% for i=1:10:length(time)
%     plot(dE_minor(:,i),z)
%     hold on
% end
% ylim([0,821])
% ylabel('ice thickness')
% xlabel('dE')

%%

C=makeConstants();


Rbottom = 1;
Rtop = 1;


H = 1000;
hw = 800;
E = 5; %creep enhancement factor
dz = 1; 
dt = 300; %(s)
z = (0:dz:H)';
wet = locatewater(hw,z);
include_ice_temperature = true;
fR_wet_variable = false;
fR_oc_variable = false;
relative_roughness_wet = NaN;
relative_roughness_OC = NaN;
fR_wet_fixed = 0.1;
fR_oc_fixed = 0.5;
Qin = 2;
Qout = Qin;

Mr = linspace(Rbottom,Rtop,length(z))';

Ms   = 0.5 * (pi .* Mr .*Mr) + 0.5 * (pi .* Mr.^2); %moulin cross-section area
Mp   = eggperimeter(Mr, Mr);   
Dh   = (4.*(pi .* Mr .* Mr)) ./ Mp; %hydrualic diameter
Rh   = (pi.* Mr .* Mr) ./ Mp; % hydraulic radius

stress.cryo = -C.rhoi*C.g*(H-z);
stress.hydro = C.rhow*C.g*(hw-z);
stress.hydro(~wet) = 0;

stress.sigx =  0e3;
stress.sigy = 50e3;
stress.tauxy = -50e3;


T_ice = C.T0*ones(size(z));



dE = elastic(Mr,stress,C);
dC = creep(Mr,z,H,stress,T_ice,dt,E,C);

[dM, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);
[dOC, Vadd_oc ] = openchannel(hw, Qin, Mr, Mr, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_oc_variable,  relative_roughness_OC,  fR_oc_fixed);



figure(333)
% plot(dC,z)
% hold on 
% plot(dE,z)
% plot(dM,z)
% hold on
% plot(dOC,z)
plot(AA,z)


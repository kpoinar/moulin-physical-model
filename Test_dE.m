

% TEST for a single timestep

H = 1000;
hw = 800;
Rtop = 1;
Rbottom = 1;

dz = 1;
z = (0:dz:H)';
Mr = linspace(Rbottom,Rtop,length(z))';


%% dE


C=makeConstants();

wet = locatewater(hw,z);

stress.cryo = -C.rhoi*C.g*(H-z);
stress.hydro = C.rhow*C.g*(hw-z);
stress.hydro(~wet) = 0;

stress.sigx =  0e3;
stress.sigy = 50e3;
stress.tauxy = -50e3;

dE = elastic(Mr,stress,C);
figure(1)
plot(dE,z)

%% dC


C=makeConstants();


E = 5;
dt = 300;

wet = locatewater(hw,z);
T_ice = C.T0*ones(size(z));

stress.cryo = -C.rhoi*C.g*(H-z);
stress.hydro = C.rhow*C.g*(hw-z);
stress.hydro(~wet) = 0;

dC = creep(Mr,z,H,stress,T_ice,dt,E,C);
figure(2)
plot(dC,z)


%% dM


dt = 300;
Qout = 2;

wet = locatewater(hw,z);
T_ice = C.T0*ones(size(z));

Ms   = 0.5 * (pi .* Mr .*Mr) + 0.5 * (pi .* Mr.^2); %moulin cross-section area
Mp   = eggperimeter(Mr, Mr);   
Dh   = (4.*(pi .* Mr .* Mr)) ./ Mp; %hydrualic diameter
Rh   = (pi.* Mr .* Mr) ./ Mp; % hydraulic radius

include_ice_temperature = true;
fR_wet_variable = false;
relative_roughness_wet = NaN;
fR_wet_fixed = 0.1;

[dM, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);
figure(3)
plot(dM,z)

%% dOC


dt = 300;
Qin = 2;
wet = locatewater(hw,z);
T_ice = C.T0*ones(size(z));

include_ice_temperature = true;
fR_oc_variable = false;
relative_roughness_OC = NaN;
fR_oc_fixed = 0.5;

[dOC, Vadd_oc ] = openchannel(hw, Qin, Mr, Mr, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_oc_variable,  relative_roughness_OC,  fR_oc_fixed);
figure(4)
plot(dOC,z)


%% dF simple


C=makeConstants();

dt = 300;
t = 300;
T_ice = 260*ones(size(z));


dF=refreeze_simple(Mr,Mr,T_ice,z,hw,t,dt,C);
figure(5)
plot(dF,z)

%% loop

C=makeConstants();


E = 5; %creep enhancement factor
 
dt = 300; %(s)
t = 9000;

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


T_ice = 260*ones(size(z));



dE = elastic(Mr,stress,C);
dC = creep(Mr,z,H,stress,T_ice,dt,E,C);

[dM, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);
[dOC, Vadd_oc ] = openchannel(hw, Qin, Mr, Mr, Mr, dt, T_ice, dz, z, wet, include_ice_temperature, fR_oc_variable,  relative_roughness_OC,  fR_oc_fixed);

dF_simple=refreeze_simple(Mr,Mr,T_ice,z,hw,t,dt,C);

figure(333)
% plot(dC,z)
% hold on 
% plot(dE,z)
% plot(dM,z)
% hold on
% plot(dOC,z)
plot(dF_simple,z)


C=makeConstants();
H = 500; %m
L = 15e3; %m
alpha = 0.0;
Rtop = 3;
Rbottom = 3;
E = 5; %creep enhancement factor
stress.sigx =  0e3;
stress.sigy = 50e3;
stress.tauxy = -50e3;

dz = 1;
z = (0:dz:H)';

dt =300;
tmax = 10*3600*24;
time = 0:dt:tmax;
% 
% Qin_mean = 1;
% dQ = 0.2;
% period = 3600*24;
% Qin = dQ * sin(2*pi*time/period) + Qin_mean;
% Qin = [0,0,Qin];
Q_pira_18 = readtable('smooth_melt_data_celia/lowc18_Q_pira_rolling.csv');
Qin = interp1(Q_pira_18.Seconds,Q_pira_18.melt_rate,time);


Ttop = C.T0;
Tbottom = C.T0;
T_ice = linspace(Tbottom,Ttop,length(z))';

include_ice_temperature = true;
fR_wet_variable = false;
fR_oc_variable = false;
relative_roughness_wet = NaN;
relative_roughness_OC = NaN;

fraction_pd_melting = 0;
fR_wet_fixed = 0.1;
fR_oc_fixed = 0.5;


%INITIALIZE

cc = 1;
%hw=zeros(length(time));
hw(1) = H; 
S = 0.5;
Qadd = 0;
Mr_major = linspace(Rbottom,Rtop,length(z))';
Mr_minor = linspace(Rbottom,Rtop,length(z))';
Mx_upstream = -Mr_major;
Mx_downstream = Mr_minor;

%LOOP
for t=time(1:end-1)
    cc = cc+1;
    
    %update moulin geom
    Ms   = 0.5 * (pi .* Mr_minor .*Mr_major) + 0.5 * (pi .* Mr_minor.^2); %moulin cross-section area
    Mp   = eggperimeter(Mr_minor, Mr_major);   
    Dh   = (4.*(pi .* Mr_minor .* Mr_major)) ./ Mp; %hydrualic diameter
    Rh   = (pi.* Mr_minor .* Mr_major) ./ Mp; % hydraulic radius
    
    %calculate head
    tspan = [0,dt];
    y0    = [hw(cc-1), S];
    opt   = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));
    Qin_compensated = Qin(cc);%+Qadd;
    [hw(cc),S,Qout, dydt_out]   = subglacialsc(Ms,z,Qin_compensated,H,L,C,dt,tspan,y0, opt);
    
    
    %update variables
    wet= locatewater(hw(cc),z);
    stress.cryo = -C.rhoi*C.g*(H-z);
    stress.hydro = C.rhow*C.g*(hw(cc)-z);
    stress.hydro(~wet) = 0;
    
    %calculate deltas
    dC_minor = creep(Mr_minor,z,H,stress,T_ice,dt,E,C);
    dC_major = creep(Mr_major,z,H,stress,T_ice,dt,E,C);

    dE_minor = elastic(Mr_minor,stress,C);
    dE_major = elastic(Mr_major,stress,C);
    
    [dM, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, Mx_downstream, dt, T_ice, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);        
    [dOC, Vadd_oc ] = openchannel(hw(cc), Qin(cc), Mr_minor, Mr_major, Mx_upstream, dt, T_ice, dz, z, wet, include_ice_temperature, fR_oc_variable,  relative_roughness_OC,  fR_oc_fixed);
    dP = potentialdrop(Qin(cc),wet,Mp,dt,C,fraction_pd_melting);
    
    dG = deformGlen(H, T_ice, alpha, z, C.n, dt, C);
    
   %Vadd_C = calculate_dQ_deformation(dC_major,dC_minor,M,z,wet);
   %Vadd_E = calculate_dQ_deformation(dE_major,dE_minor,M,z,wet);
    Qadd_tot = Vadd_turb+Vadd_oc;
    
    dr_major = dC_major+dE_major+dM;%+dOC;
    dr_minor = dC_minor+dE_minor+dM;
    
    %update moulin shape
    Mr_major = Mr_major+dr_major;
    Mr_minor = Mr_minor+dr_minor;
    
    Mx_upstream = Mx_upstream + dG - dr_major;
    Mx_downstream = Mx_downstream + dG +dr_minor;
    
    
    
    
end
%%
figure(1)
plot(Mx_upstream,z)
hold on
plot(Mx_downstream,z)
xlim([-3,3])

figure(2)
plot(time,hw)




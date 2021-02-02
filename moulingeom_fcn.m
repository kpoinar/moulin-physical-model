% Solve for moulin geometry, considering
function time = moulingeom_fcn( workingdirectory, savelocation, makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savefile) 

% 1. elastic deformation (opening / closure)
% 2. creep deformation (opening / closure)
% 3. refreezing (closure)
% 4. turbulent melting (opening)
% 5. open-channel melting (opening)
% 6. xz-shear deformation
%
%


    cd(workingdirectory)
    
    % Do you want to plot model results?
    make_simple_plots =  makeplots_tf;
    save_figures      =  savefigures_tf;
    visible_figures   =  showfigures_tf;

    % What type of open channel do you want to do?
    %1 --> open channel parameterization
    %2 --> waterfall like parameterization
    %3 --> potential drop parameterization
    %=/1-3 --> no melting above the water line
    unfilled_melting = modelinputs.unfilled_melting; %1;
    
    %% define some basic parameters
    C         = makeConstants;  %constants used for parameterizations  %if there is a particular value we want to change gradually, 
    if isfield(modelinputs,'A_value')%, 'var')
        C.A = modelinputs.A_value; 
        C.c2 = 1*C.A*C.n^(-C.n); % Need to recalculate the closure parameter (Schoof 2010) with the new A
    end
    
    if isfield(modelinputs,'ShearModulus')%, 'var')
        C.Y = modelinputs.ShearModulus; 
    end
    
    if isfield(modelinputs,'f_sub')%, 'var')
        C.f = modelinputs.f_sub;
    end
    
    Tdatatype = modelinputs.Tdatatype{1}; %'Ryser_foxx';   %ice temperature profile to extrapolate from
    numofdays = modelinputs.numofdays; % 3  %set the number of days for the model run
    H         = modelinputs.H;         % 800   % ice thickness, meters
    R0        = modelinputs.R0;        % 2 m   % radius @bed of moulin initially
    Rtop      = modelinputs.Rtop;      % 0.5 m % radius @sfc of moulin initially
    L         = modelinputs.L;         % L  % Length of the subglacial channel
    f         = modelinputs.fract_pd_melting;         %f   % fraction of the potential energy used to open the top of the moulin (above water level)
    alpha     = modelinputs.alpha;     % 0.03      % regional surface slope (unitless), for use in Glen's Flow Law
    n         = 3;              % flow law exponent (Glen's Flow Law)
    
    %inital guesses for subglacial model
    hw(1) = H;                  % moulin water level (m)
    S(1)  = R0;%1.5* R0;                 % subglacial channel cross sectional area (m^2)
    chebx     = 0;              % chebx=1 is not working yet
    E         = modelinputs.E; %5;             % enhancement factor for creep
    
    %% set the vertical model components
    dz        = 1; %  vertical spacing, meters
    z         = (0:dz:H)';
    time.z    = z; %save the z profile in time
    %% set the duration of the model run
    sec       = 86400*numofdays;   %seconds * days
    dt        = 300;        % Timestep, seconds 
    tmax      = sec;        % seconds (5 years here)
    time.t    = dt:dt:tmax; % seconds
    time.dt   = dt;
    time.sec  = sec;
    %% set Qin
    % Construct an approximate Qin for each timestep, based on air temps:
    % Qin = double(Tair>C.To) .* (Tair-C.To)*const;
    % Qin = zeros(size(time.t));
    % Use predetermined Qins of various types
    qinreal = modelinputs.Qinreal;
    
    if qinreal 
        load(modelinputs.Qinfile{1}) %1 = time, 2 cosine function, 3
        Q(:,1) = Q2.(modelinputs.Qin_year{1}).time_seconds;
        Q(:,2) = modelinputs.Qin_dampen .* Q2.(modelinputs.Qin_year{1}).(modelinputs.Qin_basin{1});
        Qin     = interp1(Q(:,1), Q(:,2), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
        
        %apply a smoothing to Qin to dampen diurnal varibility
        Qin = max(0.1, Qin);
        Qin = smoothdata(Qin, 'movmean',  modelinputs.Qin_smoothval);
        time.Qin = Qin;  %save for future plotting
        clear Q
        
        % apply a base flow only to the subglacial system - this kind of
        % mimics a number of processes
        load(modelinputs.Qinfile{1});
        baseflow  =  Q2.(modelinputs.Qin_year{1}).(modelinputs.Qin_baseflow{1});

        Qbase  = 4.* interp1(baseflow(:,1), baseflow(:,2), time.t, 'spline', 'extrap'); 
        Qbase = Qbase';
        time.Qbase = Qbase;
        clear Qbase_subglacial baseflow

%lca commented these out in order to prevent Qbase from contributing to moulin geometry          
%         Qin = Qin + Qbase;
%         time.Qin = Qin + Qbase;
        
    else 
        load(modelinputs.Qinfile{1}) %1 = time, 2 cosine function, 3
        Qin     = interp1(Q(:,1), Q(:,modelinputs.Qintype), time.t, 'spline', 'extrap'); % run an interp just in case the timeframe changes
        %not sure if we should really have this... might require an adjustment
        %in the qin values
        %Qin     = Qin*modelinputs.Qinmultiplier1 +modelinputs.Qinmultiplier2; %scale Qin to deal with a few model issues
        % Normalize Qin
        Qin = (Qin - mean(Qin)) / range(Qin);
        % Rescale using user-defined Qin parameters
        Qin = max(0.1, Qin * modelinputs.Qinrange + modelinputs.Qinbase);
        time.Qin = Qin;  %save for future plotting
        clear Q 
        
        Qbase = time.t .* 0 ; %apply zero baseflow to the subglacial model
        
    end
    
figure
hold on 

plot(time.t, Qin)
%yyaxis right 


    
    %% set Ice temperature characteristics
    
    Tfar    = importTz(Tdatatype,z); % Kelvin
    xmax    = 30;% 80; % meters; how far away from moulin to use as infinity
    [x,dx,nx] = setupx(dt,chebx,xmax,C);
    T       = Tfar*ones(size(x));  % Ambient ice temperature everywhere to start
    T(:,1)  = C.T0;   % Melting point at the moulin wall
    time.icetemp = Tfar; %just save in the time file for reference
    %% define initial moulin characteristics
    Mrmin   = 1e-9;  % 1 mm
%     M.r_minor = R0*ones(size(z)); %To use this, the moulin should be filled
%     M.r_major = R0*ones(size(z)); %To use this, the moulin should be filled
    M.r_minor = linspace(R0,Rtop,length(z))'; % Start as a "tree" moulin with 
    M.r_major = M.r_minor;                   % R= R0 at the bottom and R= 0.1 m at the top
                                        % The geom at top can only
                                        % grow (open channel), never
                                        % shrink (no stresses).
    
    % initalize the horizontal coordinate system
    %This assumes that ice flow is from left to right
    M.xu = -M.r_major;
    M.xd =  M.r_minor;
    % Pin the bed of the upstream wall to x=0 while retaining the initial
    % moulin shape / radius:
    x0 = M.xu(1);
    M.xu = M.xu - x0;
    M.xd = M.xd - x0;
    %initalize dVdt
    Qadd = 0;
    
    %% Set turbulence parameters
   
    fR_wet_variable = modelinputs.variable_fR_wet;
    relative_roughness_wet = modelinputs.relative_roughness; %0.2; %increasing this value increases the amount of melting due to turbulence.
    fR_wet_fixed       = modelinputs.fR_fixed_wet;
    
    fR_oc_variable = modelinputs.variable_fR_oc;
    relative_roughness_OC = modelinputs.relative_roughness_OC; %1e-9;%1e-12;  % This one modifies the melt from open channel flow.
    fR_oc_fixed           = modelinputs.fR_fixed_oc;
    
 
    
    include_ice_temperature = true; %true means that the change in the ice temperature is included in...
    %the calculated change in moulin radius. If false, it makes the implicit
    %assumption that the ice temperature and water temperature are both at the pressure melting temperature.
    
    if include_ice_temperature
        Ti = Tfar;
    else
        Ti = NaN; %#ok<UNRCH>
    end
    
     %true means that the friction factor is calculated using..
    %the Bathurst formulation... this equation is valid when
    %roughness height ./ hydrualic diameter >= 0.05
    % if false, the Colebrook-White formulation will be applied, which is only
    % valid when roughness height ./ hydrualic diameter < 0.05
    
    %% Assign elastic deformation parameters
    stress.sigx = modelinputs.sigx;
    stress.sigy = modelinputs.sigy;
    stress.tauxy = modelinputs.tauxy;
    
    %% save general parameters in time file
    time.parameters.stress = stress;
    time.parameters.relative_roughness = relative_roughness_wet;
    time.parameters.relative_roughness_OC = relative_roughness_OC;
    time.parameters.creepenhancement = E;
    time.parameters.H = H;
    time.parameters.L =L;
    time.parameters.R0 = R0;
    time.parameters.numofdays =  numofdays;
    time.parameters.f = f;
    
    if unfilled_melting ==1
        time.parameters.dOCtype = 'open channel melting parameterization';
    elseif unfilled_melting ==2
        time.parameters.dOCtype = 'waterfall like melting parameterization';
    elseif unfilled_melting ==3
        time.parameters.dOCtype = 'potential drop melting parameterization';
    else
        time.parameters.dOCtype = 'no melting applied above waterline';
    end
    
    
    %% Step through time
    cc = 0;
    nt = length(time.t);
    for t = time.t


        
        
        
        
        cc = cc+1;
        if ~mod(cc,100), fprintf('timestep %d of %d (%1.0f%%) after %1.1f minutes \n', cc,nt,cc/nt*100,toc/60); end
        % Consider using the previous moulin radius in all calculations in each
        % timestep, so that the final result is not dependent on the order in
        % which I do creep, refreeze, turbulent melt, elastic, etc.
        Mrminor_prev  = M.r_minor;
        Mrmajor_prev  = M.r_major;
        Mxuprev = M.xu;

        

        
        
        %%%%%%%%%%
        %calculate moulin parameters
        Ms   = 0.5 * (pi .* Mrminor_prev .*Mrmajor_prev) + 0.5 * (pi .* Mrminor_prev.^2); %moulin cross-section area
        Mp   = eggperimeter(Mrminor_prev, Mrmajor_prev);   %pi.* (3 .*(Mrminor_prev + Mrmajor_prev) - sqrt((3.* Mrminor_prev + Mrmajor_prev) .* (Mrminor_prev +3 .* Mrmajor_prev))); % wetted/melting perimeter =  ellipse perimeter approx pi [ 3(Mrminor+Mrmajor) - sqrt((3*Mrminor+Mrmajor)(Mrminor+3*Mrmajor))]
        Dh   = (4.*(pi .* Mrminor_prev .* Mrmajor_prev)) ./ Mp; %hydrualic diameter
        Rh   = (pi.* Mrminor_prev .* Mrmajor_prev) ./ Mp; % hydraulic radius
        
        Ms_prev = Ms;
        
        %%%%%%%%%%
        % Subglacial Schoof model: Conduit size
        tspan = [0,dt];
        y0    = [hw, S];
        %[hw,S,Qout]   = subglacialsc(Mrminor_prev,z,Qin(cc),H,L,C,tspan,y0);
        opt   = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));
        %Qin_tot       = Qin(cc) + time.V
        
        %Qin_compensated = Qin(cc)+Qadd;
        Qin_compensated = Qin(cc)+Qadd + Qbase(cc); %including Qbase provides a minimum base flow for the subglacial model without needing to route additional water through the moulin
        time.Qin_compensated(cc) = Qin_compensated;
        %[hw,S,Qout]   = subglacialsc(Mrminor_prev,z, Qin_subbase(cc),H,L,C,tspan,y0, opt); %consider adding Vadd to the qin values
        [hw,S,Qout, dydt_out]   = subglacialsc(Ms,z,Qin_compensated,H,L,C,dt,tspan,y0, opt); %consider adding Vadd to the qin values
            %the first term in the function had been MrMinor_prev, which
            %actually is the radius, not the cross-sectional area, LCA
            %fixed on 6/7/20

            time.S(cc)    = S;
            time.hw(cc)   = hw;
            time.Qout(cc) = Qout ;% - Qbase(cc); %this removes the baseflow from the Qout
            time.dydt_out(cc) = dydt_out;
            time.Qbase(cc) = Qbase(cc); 
        
        %%%%%%%%%%%
        % which nodes are underwater or at the water line (wet) versus above the water line?
        wet = locatewater(hw,z);
        time.wet(:,cc)    = wet;
        % Calculate hydrostatic pressures everywhere
        % Ice hydrostatic stress (INWARD: Negative)
        stress.cryo = -C.rhoi*C.g*(H-z);
        % Water hydrostatic stress (OUTWARD: Positive)
        stress.hydro = C.rhow*C.g*(hw-z);
        stress.hydro(~wet) = 0; % Anywhere that is not wet does not have the opening force from water
        
        
        %%%%%%%%% dC: Creep deformation
        %Creep deformation: do this first because it is a larger term
        dC_minor = creep(Mrminor_prev,z,H,stress,T,dt,E,C);
            time.dC_minor(:,cc) = dC_minor;
        dC_major = creep(Mrmajor_prev,z,H,stress,T,dt,E,C);
            time.dC_major(:,cc) = dC_major;
        % Creep deformation changes the volume of the moulin.  We need to 
        % calculate this volume change and send it to subglacialsc:
        Vadd_C = calculate_dQ_deformation(dC_major,dC_minor,M,z,wet);
        %time.Vadd_C(cc) = Vadd_C;
        
        %%%%%%%%% dF: Refreezing
        % Refreezing
        %     T(z>hw,1) = Tair(cc);
        %     [~,dF,T,Vfrz] = refreeze(Mrminor_prev,T,z,hw,wet,dF,nx,x,dx,dt,C);
        % 	  [~,dF,Vfrz]   = refreeze_simple(Mrminor_prev,z,)
        
        
        
        %%%%%%%%% dM: Turbulent melting
        % Turbulent melting:
        [dM, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, M.xd, dt, Ti, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);
        time.dM(:,cc)  =  dM;
        time.uw(:,cc)  =  uw;
        % time.V(cc)  = Vadd;
        %[dM_major, uw_major] = turbulence(hw, Qout, Mrmajor_prev, dt, Ti, z, relative_roughness, Bathurst, include_ice_temperature);
        %    time.dM_major(:,cc)  =  dM_major;
        %    time.uw_major(:,cc)  =  uw_major;
        %time.Vadd_major(cc)  = Vadd_major;
        % Calculate the water volume added to the moulin by calculating the enlargement of the moulin due to turbulent melting
        %Vadd_turb = waterVolumeFromTurbulence(Mrminor_prev, Mrmajor_prev, dM, z, wet);
        %        Vadd_turb = waterVolumeFromTurbulence(Mrminor_prev, Mrmajor_prev, dM, z, wet, dt);
            time.Vadd_turb(cc) = Vadd_turb;
        
        
        %%%%%%%%% dOC: Melting due to open channel flow above the moulin water line
            [dOC, Vadd_oc ] = openchannel(hw, Qin(cc), Mrminor_prev, Mrmajor_prev, M.xu, dt, Ti, dz, z, wet, include_ice_temperature, fR_oc_variable,  relative_roughness_OC,  fR_oc_fixed);
            
            % Scale the open channel displacement down by 1/2 to reflect the
            % displacement at exactly the upstream point:
            %dOC = dOC / 2; %lca commented this out after discussion with
            %KP on 4/27.
            
            time.dOC(:,cc)  =  dOC;
            time.Vadd_oc(cc)    =  Vadd_oc;
%             time.oc_dL(:,cc) = dL;
%             time.oc_hL(:,cc) = hL;
%             time.oc_fR(:,cc) = fR;
%             time.oc_neg(:,cc) = neg;
%             time.oc_rh(:,cc) = rh;
%             time.oc_cres_area(:,cc) = cres_area;
%             time.oc_Mp(:,cc) = Mp;
        %%%%%%%%% dP: Expansion from gravitational potential energy above the water
            %%%%%%%%% line
            dP = potentialdrop(Qin(cc),wet,Mp,dt,C,f);
            % The reason for calculating the above is to offset the elastic closure
            % at the top of the moulin.  On its own, the moulin will close
            % elastically after some days to months (depending on C.Y).  We know
            % that does not happen.  Hence, we add some turbulent "waterfall"
            % melting above the water line.
            time.dP(:,cc) = dP;
            time.Vadd_p   = 0;
            Vadd_p = 0;       
        
        
        %%%%%%%%% Vadd: Added water from melted ice into Qin
        % NOTE 2 MARCH 2020: This is a large amount of meltwater (~8 m2 per dt)
        % compared to the current Qin (~4 m2 per dt).  It can break the model
        % if the moulin and subglacial conduit aren't big enough.
        % Add the
        %CT commented below: Qin is update just before subglacialsc
        %if cc < length(time.t)
        %    Qin(cc+1) = Qin(cc+1) + Vadd_turb / dt + Vadd_oc / dt + Vadd_p / dt;
        %end
        
        %%%%%%%%% dE: Elastic deformation
        % Elastic deformation: This is small, and sensitive to water pressure
        dE_minor = elastic(Mrminor_prev,stress,C);
        time.dE_minor(:,cc) = dE_minor;
        dE_major = elastic(Mrmajor_prev,stress,C);
        time.dE_major(:,cc) = dE_major;
        % Elastic deformation changes the volume of the moulin.  We need to 
        % calculate this volume change and send it to subglacialsc:
        Vadd_E = calculate_dQ_deformation(dE_major,dE_minor,M,z,wet);
        %time.Vadd_E(cc) = Vadd_E;
        
        %%%%%%%%% dG: Asymmetric deformation due to Glen's Flow Law
        dG = deformGlen(H, T, alpha, z, n, dt, C);
        time.dG(:,cc) = dG;
        
        
        % Update Qadd, the additional "water flux" from elastic/creep
        % deformation (this is phantom water) and open-channel flow above
        % the water line (this is real water):
        Qadd=(Vadd_E+Vadd_C+Vadd_oc)/dt; % Divide volumes by timestep to get Q
        time.Qadd(cc) = Qadd;
        
        
        %%%%%%%%LCA March 24 --> these need to change to reflect new dOC so
        %%%%%%%%both the  xd and xu need to have dOC depending on the dOC
        %%%%%%%%choices
        %%%%%%%%parameterizations
        % Calculate the horizontal position of the moulin within the ice column
        
        
        %This if statement allows you to determine if you want potential
        %drop applied on the downstream radius
        if modelinputs.use_pD_downstream
            M.xu = M.xu - dC_major - dE_major - dM + dG - dOC - 0 * dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
            % Important Note: the +dG above is correct.
            % The upstream wall moves downstream.
            
            M.xd = M.xd + dC_minor + dE_minor + dM + dG + 0 * dOC +  dP; % if you dont want any melt from dP, then use 0 * dP
            % The downstream wall also moves downstream
            % at the same rate, dG.
            
            % Shift them both back upstream so that the bed of the upstream wall
            % stays pinned at x = 0:
            x0 = M.xu(1);
            M.xu = M.xu - x0;
            M.xd = M.xd - x0;
            %
            % Now use the moulin positions to assign the major and minor radii:
            M.r_minor = max(M.r_minor + dC_minor + dE_minor + dM + dP, Mrmin); %not sure whhat the max does here, but added dP
            M.r_major = (M.xd - M.xu) - M.r_minor;
            
            
        else
            
            M.xu = M.xu - dC_major - dE_major - dM + dG - dOC - 0*dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
            % Important Note: the +dG above is correct.
            % The upstream wall moves downstream.
            
            M.xd = M.xd + dC_minor + dE_minor + dM + dG + 0*dOC +  0*dP; % if you dont want any melt from dP, then use 0 * dP
            % The downstream wall also moves downstream
            % at the same rate, dG.
            
            % Shift them both back upstream so that the bed of the upstream wall
            % stays pinned at x = 0:
            x0 = M.xu(1);
            M.xu = M.xu - x0;
            M.xd = M.xd - x0;
            %
            % Now use the moulin positions to assign the major and minor radii:
            M.r_minor = max(M.r_minor + dC_minor + dE_minor + dM + 0 * dP, Mrmin); %not sure whhat the max does here, but added dP
            M.r_major = (M.xd - M.xu) - M.r_minor;
            

            
        end
 
        

        

        
        % Record the used moulin geometry
        time.M.r_minor(:,cc) = M.r_minor;
        time.M.r_major(:,cc) = M.r_major;
        time.M.xu(:,cc) = M.xu;
        time.M.xd(:,cc) = M.xd;
        %
        % Record volume capacity of moulin
        % This reflects the semi-ellipse, semi-circular geometry
        % and uses the variable "wet" for where is water.
        [time.Mcapacity(cc), time.Wvolume(cc)] = moulincapacity(M,z,wet);



        
    end
    %% figures
    
     if make_simple_plots
         simpleplots(time, save_figures, savelocation, visible_figures, modelinputs.runnumber, savefile)
     end
    
    %%
    
 







end

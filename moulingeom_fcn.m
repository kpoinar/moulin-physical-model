% Solve for moulin geometry, considering
function [time, elasticcomp] = moulingeom_fcn( workingdirectory, savelocation, makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savefile) 

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
    
    %Define constants, including depedents
    C         = makeConstants;  %constants used for parameterizations  %if there is a particular value we want to change gradually, 
    
 
    if isfield(modelinputs,'A_value')%, 'var')
        C.A = modelinputs.A_value; 
        C.c2 = 2*C.A*C.n^(-C.n); % Need to recalculate the closure parameter (Schoof 2010) with the new A
    end
    
    if isfield(modelinputs,'ShearModulus')
        C.Y = modelinputs.ShearModulus; 
    end
    
    if isfield(modelinputs,'f_sub')
        C.f = modelinputs.f_sub;
    end
    
    Tdatatype = modelinputs.Tdatatype{1};    % ice temperature profile to extrapolate from
    numofdays = modelinputs.numofdays;       % set the number of days for the model run
    H         = modelinputs.H;               % ice thickness, meters
    R0        = modelinputs.R0;              % radius @bed of moulin initially
    Rtop      = modelinputs.Rtop;            % radius @sfc of moulin initially
    L         = modelinputs.L;               % Length of the subglacial channel
    f         = modelinputs.fract_pd_melting;% fraction of the potential energy used to open the top of the moulin (above water level)
    alpha     = modelinputs.alpha;           % regional surface slope (unitless), for use in Glen's Flow Law
    n         = 3;                           % flow law exponent (Glen's Flow Law)
    
    %inital guesses for subglacial model
    hw(1) =  H;                              % moulin water level (m)
    S(1)  = R0;                              % subglacial channel cross sectional area (m^2)
    chebx     = 0;                           % chebx=0 is the only option (0,1) that has been throughly tested
    E         = modelinputs.E;               % enhancement factor for creep
    
    %% set the vertical model components
    dz        = 1;                          % vertical spacing, meters
    z         = (0:dz:H)';                  % vertical grid 

    %% set the duration of the model run
    sec       = 86400*numofdays;            % seconds * days
    dt        = 300;                        % Timestep, seconds (longer timesteps work but can lengthen the duration of the run due to ODE stiffness)
    tmax      = sec;                        % maximum timestep
    time.t    = dt:dt:tmax;                 % time vector in seconds saved in outputs


    %% set Qin
    %Qin is the discharge entering the moulin
    %Qbase is the additional amount of water added to the subglacial system
           %in this formulation, Qbase is a smoothed function of Qin. It is
           %only applied during the realistic scenarios and is calculated
           %as:         Qbase  = baseflowmultiplier.* interp1(baseflow(:,1), baseflow(:,2), timet, 'spline', 'extrap'); 
    
    baseflowmultiplier = 2; 
   
    [Qin, Qbase ] = Qincalc(modelinputs.Qinreal, modelinputs.Qinfile,...
                    modelinputs.Qin_smoothval, modelinputs.Qin_year, modelinputs.Qin_basin, modelinputs.Qin_baseflow,...
                    modelinputs.Qin_dampen, modelinputs.Qintype, modelinputs.Qinrange,...
                    modelinputs.Qinbase, time.t, baseflowmultiplier);
                
                

    
    %% set Ice temperature characteristics
    
    Tfar    = importTz(Tdatatype,z);    % Kelvin
    xmax    = 30;                       % meters; how far away from moulin to use as infinity
    [x,~,~] = setupx(dt,chebx,xmax,C);
    T       = Tfar*ones(size(x));       % Ambient ice temperature everywhere to start
    T(:,1)  = C.T0;                     % Melting point at the moulin wall


    
    %% define initial moulin characteristics
    Mrmin   = 1e-9;  % 1 mm

    %Initialize the vertical coordinate system.
    M.r_minor = linspace(R0,Rtop,length(z))'; % This sets the initial moulin radius. 
    M.r_major = M.r_minor;                    % if R0 and Rtop are different than 
                                              %   the difference is linearly interpolated.
                                              % The geom at top can only
                                              % grow (open channel), never
                                              % shrink (no stresses).
    
    %Initalize the horizontal coordinate system. 
    %This assumes that ice flow is from left to right
    M.xu = -M.r_major;
    M.xd =  M.r_minor;
    
    % Pin the bed of the upstream wall to x=0 while retaining the initial
    % moulin shape / radius:
    x0 = M.xu(1);
    M.xu = M.xu - x0;
    M.xd = M.xd - x0;
    
    %Initalize added discharge associated with melting. (dV/dt)
    Qadd = 0;
    
    %% Set turbulence parameters
   
    % Initialize the roughness for water filled portion of the moulin    
    relative_roughness_wet  = modelinputs.relative_roughness; %increasing this value increases the amount of melting due to turbulence.
    fR_wet_variable         = modelinputs.variable_fR_wet;
    fR_wet_fixed            = modelinputs.fR_fixed_wet;
    
    % Initialize the roughness of the non-water filled portion of the
    % moulin
    relative_roughness_OC   = modelinputs.relative_roughness_OC; % This one modifies the melt from open channel flow.
    fR_oc_variable          = modelinputs.variable_fR_oc;
    fR_oc_fixed             = modelinputs.fR_fixed_oc;
    
    % Include ice temperature in the calculation of wall melting or not.
    %true means that the change in the ice temperature is included in...
    %the calculated change in moulin radius. If false, it makes the implicit
    %assumption that the ice temperature and water temperature are both at the pressure melting temperature.
    include_ice_temperature = true; 
    
    if include_ice_temperature
        Ti = Tfar;
    else
        Ti = NaN; %#ok<UNRCH>
    end
    
    %% Assign elastic deformation parameters
  
    stress.sigx = modelinputs.sigx;
    stress.sigy = modelinputs.sigy;
    stress.tauxy = modelinputs.tauxy;
    
    %% save general parameters in time (output) file
    
    time.parameters.stress                = stress;
    time.parameters.relative_roughness    = relative_roughness_wet;
    time.parameters.relative_roughness_OC = relative_roughness_OC;
    time.parameters.creepenhancement      = E;
    time.parameters.H                     = H;
    time.parameters.L                     = L;
    time.parameters.R0                    = R0;
    time.parameters.numofdays             =  numofdays;
    time.parameters.f = f;
    time.icetemp = Tfar;                %just save in the time file for reference
    time.Qin = Qin;
    time.Qbase = Qbase;  
    time.dt   = dt;                         % dt saved in outputs
    time.sec  = sec;  
    time.z    = z;                          % save the z profile in time structure.
          
    %This statement is designed to save clear information on what type of melting is prescribed above the water line.      
    if unfilled_melting ==1
        time.parameters.dOCtype = 'open channel melting parameterization';
    elseif unfilled_melting ==2
        time.parameters.dOCtype = 'waterfall like melting parameterization';
    elseif unfilled_melting ==3
        time.parameters.dOCtype = 'potential drop melting parameterization';
    else
        time.parameters.dOCtype = 'no melting applied above waterline';
    end
    
    
    %% initialize variable for elastic2 
        dR_major_total = zeros(length(M.r_major), 1);
        dR_minor_total = zeros(length(M.r_major), 1);
    
    
    %% Step through time
    cc = 0;
    nt = length(time.t);
    for t = time.t
        
        cc = cc+1;
        if ~mod(cc,100), fprintf('timestep %d of %d (%1.0f%%) after %1.1f minutes \n', cc,nt,cc/nt*100,toc/60); end
        % Use the previous moulin radius in all calculations in each
        % timestep, so that the final result is not dependent on the order
        % of creep, refreeze, turbulent melt, elastic, etc. calculations
        Mrminor_prev  = M.r_minor;
        Mrmajor_prev  = M.r_major;
        Mxuprev       = M.xu;
        hw_prev       = hw;
        
        
        %%%%%%%%%%
        %calculate moulin hydraulic parameters based on the chosen moulin
        %plan shape
        if contains(modelinputs.planshape, 'egg')
            Ms   = 0.5 * (pi .* Mrminor_prev .*Mrmajor_prev) + 0.5 * (pi .* Mrminor_prev.^2); %moulin cross-section area
            Mp   = eggperimeter(Mrminor_prev, Mrmajor_prev);   %pi.* (3 .*(Mrminor_prev + Mrmajor_prev) - sqrt((3.* Mrminor_prev + Mrmajor_prev) .* (Mrminor_prev +3 .* Mrmajor_prev))); % wetted/melting perimeter =  ellipse perimeter approx pi [ 3(Mrminor+Mrmajor) - sqrt((3*Mrminor+Mrmajor)(Mrminor+3*Mrmajor))]
            Dh   = (4.*(pi .* Mrminor_prev .* Mrmajor_prev)) ./ Mp; %hydrualic diameter
            Rh   = (pi.* Mrminor_prev .* Mrmajor_prev) ./ Mp; % hydraulic radius
        
            Ms_prev = Ms;
            
        elseif contains(modelinputs.planshape, 'circle') % Approximate the circle as an ellipse because this allows us to keep a major and minor radius. Need to confim
            Ms   =  pi .* Mrmajor_prev .* Mrmajor_prev; %moulin cross-section area as a circle
            Mp   = 2 .* pi .* Mrmajor_prev; %pi.* (3 .*(Mrminor_prev + Mrmajor_prev) - sqrt((3.* Mrminor_prev + Mrmajor_prev) .* (Mrminor_prev +3 .* Mrmajor_prev))); % wetted/melting perimeter =  ellipse perimeter approx pi [ 3(Mrminor+Mrmajor) - sqrt((3*Mrminor+Mrmajor)(Mrminor+3*Mrmajor))]
            Dh   = (4.*(pi .* Mrmajor_prev .* Mrmajor_prev)) ./ Mp; %hydrualic diameter
            Rh   = (pi.* Mrmajor_prev .* Mrmajor_prev) ./ Mp; % hydraulic radius
        
            Ms_prev = Ms;
             
            
        end
        
        
%%%%%%Subglacial model: Conduit size%%%%%%%%%%%%%  
        
        % Calculate the total amount of water moving through the subglcial
        % channel moulin Qin, additional melt (Qadd), and the prescribed
        % baseflow (Qbase)
        Qin_subtotal             = Qin(cc)+ Qadd + Qbase(cc); %including Qbase provides a minimum base flow for the subglacial model without needing to route additional water through the moulin
        time.Qin_subtotal(cc)    = Qin_subtotal;
        
        % Conditions needed for the ODE solver
        tspan                    = [0, dt];
        y0                       = [hw, S];
        opt                      = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3)); %first guesses for the ode solvers
        
        % ODE to solve subglacial model
        hw_init = hw;
        [hw,S,Qout, dydt_out]    = subglacialsc(Ms,z,Qin_subtotal,H,L,C,dt,tspan,y0, opt); 
        %[hw,S,Qout, dydt_out]   = subglacialsc_fixedS(Ms,z,Qin_subtotal,H,L,C,dt,tspan,y0, opt); 
        
        % Assign outputs to the output vector
        time.S(cc)               = S;
        time.hw(cc)              = hw;
        time.Qout(cc)            = Qout ; 
        time.dydt_out(cc)        = dydt_out;
        time.Qbase(cc)           = Qbase(cc); 
        
      
        
%%%%%Node pressures%%%%%%%%%%%%%
% Submerged nodes and cryostatic and hydrostatic stresses        
        
        % Determine submerged and non-submerged model nodes 
        wet               = locatewater(hw,z);
        time.wet(:,cc)    = wet;
        
        
        % Calculate hydrostatic pressures everywhere
        % Ice hydrostatic stress (INWARD: Negative)
        stress.cryo        = -C.rhoi * C.g * (H - z);
        
        % Water hydrostatic stress (OUTWARD: Positive)
        stress.hydro       = C.rhow * C.g * (hw - z);
        stress.hydro(~wet) = 0; % Anywhere that is not wet does not have the opening force from water
        
        
%%%%%dC: Creep deformation%%%%%%%%%%%%%
        %Creep deformation: do this first because it is a larger term

        dC_major                    = creep(Mrmajor_prev,z,H,stress,T,dt,E,C);
        time.dC_major(:,cc)         = dC_major;
        
        % minor radius change is dependent on  the plan shape
         if contains(modelinputs.planshape, 'egg')    
            dC_minor                = creep(Mrminor_prev,z,H,stress,T,dt,E,C);
                time.dC_minor(:,cc) = dC_minor;    
         elseif contains(modelinputs.planshape, 'circle')
            dC_minor                = dC_major;
                time.dC_minor(:,cc) = dC_minor;
         end
             
        % Creep deformation changes the volume of the moulin.  We need to 
        % calculate this volume change and send it to subglacialsc:
        Vadd_C = calculate_dQ_deformation(dC_major,dC_minor,M,z,wet);
        
        
%%%%% dF: Refreezing%%%%%%%%%%%%%
% Refreezing for outside the melt season. This is applied uniformly around
% the moulin radius
        
        % [dF, Vfrz]   = refreeze_simple(Mrminor_prev,Mrmajor_prev, Tfar,...
        %                                z, hw, time.t, dt, C);
        
        %This complex refreezing paramterization is still in development.
        %No gaurentee of functioning
        %     T(z>hw,1) = Tair(cc);
        %     [~,dF,T,Vfrz] = refreeze(Mrminor_prev,T,z,hw,wet,dF,nx,x,dx,dt,C);
        
        
        
        %%%%%%%%% dM: Turbulent melting
        % Turbulent melting:
        
        
        
        
        [dMarea, uw, Vadd_turb] = turbulence(Qout, Ms, Mp, Dh, Rh, M.xd, dt, Ti, dz, z, wet, include_ice_temperature, fR_wet_variable, relative_roughness_wet,   fR_wet_fixed);
        
        %dM is the change in cross-sectional area. Now it needs to be
        %converted into a change in radius. Under both prescribed planform
        %shapes, the change in radii is uniform across major and minor axes
        % however, the equations to derive dM 
        
        if contains(modelinputs.planshape, 'egg')
%            dM = 2 .* dMarea ./ (pi .* (5 .* Mrmajor_prev + 3 .* Mrminor_prev - sqrt((3 .* Mrmajor_prev + Mrminor_prev).*(Mrmajor_prev + 3 .* Mrminor_prev)))); 
            dM = 2 .* dMarea ./ (pi .* (5 .* Mrminor_prev + 3 .* Mrmajor_prev - sqrt((3 .* Mrminor_prev + Mrmajor_prev).*(Mrminor_prev + 3 .* Mrmajor_prev)))); 

        elseif contains(modelinputs.planshape, 'circle')
            dM = sqrt((Ms_prev + dMarea)./pi) - Mrmajor_prev;
        end
        
        
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
            
        if contains(modelinputs.planshape, 'egg')
            time.dOC(:,cc)  =  dOC;
            time.Vadd_oc(cc)    =  Vadd_oc;
        elseif contains(modelinputs.planshape, 'circle')
            %This applies the melting uniformly around the perimeter...
            time.dOC(:,cc)  =  dOC;
            time.Vadd_oc(cc)    =  Vadd_oc;
        end
            
            
            
            
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
         
  %%%%%%%%%Implement the new elastic scheme elastic2

        
         [dE_major, majordMr_dt, majordP_dt, majorMr, majorP1] = elastic2(Mrmajor_prev, stress, hw_prev, dR_major_total, dt, C,z,wet); 
         time.dE_major(:,cc) = dE_major;
        if contains(modelinputs.planshape, 'egg')
            [dE_minor, minordMr_dt, minordP_dt, minorMr, minorP1] = ...
                        elastic2(Mrminor_prev, stress, hw_prev, dR_minor_total, dt, C,z, wet); 
                time.dE_minor(:,cc) = dE_minor;
        elseif contains(modelinputs.planshape, 'circle')
                dE_minor = dE_major;
                time.dE_minor(:,cc) = dE_minor;
        end         
         
        dR_major_total = dE_major - dOC + dC_major - dM +dP;
        dR_minor_total = dE_minor - 0*dOC + dC_major - dM + dP;
        
        
        elasticcomp.dE_major(:,cc)     = dE_major;
        elasticcomp.dE_minor(:,cc)     = dE_minor;
        elasticcomp.majordMr_dt(:,cc)  = majordMr_dt;
        elasticcomp.majordP_dt(:,cc)   = majordP_dt;
        elasticcomp.majorMr(:,cc)      = majorMr;
        elasticcomp.majorP1(:,cc)      = majorP1;
        elasticcomp.minordMr_dt(:,cc)  = minordMr_dt;
        elasticcomp.minordP_dt(:,cc)   = minordP_dt;
        elasticcomp.minorMr(:,cc)      = minorMr;
        elasticcomp.minorP1(:,cc)      = minorP1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%3/22/22 commented by LCA to test elastic2        
%         dE_major = elastic(Mrmajor_prev,stress,C);
%         time.dE_major(:,cc) = dE_major;
%        
%         %This if statement is used to eliminate the minor radius when the model us run with a circular crossectional area.
%         if contains(modelinputs.planshape, 'egg')
%             dE_minor = elastic(Mrminor_prev,stress,C);
%                 time.dE_minor(:,cc) = dE_minor;
%         elseif contains(modelinputs.planshape, 'circle')
%             dE_minor = dE_major;
%                 time.dE_minor(:,cc) = dE_minor;
%         end
         
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
            M.xu = M.xu - dC_major - dE_major - dM +  dG - dOC - 0 * dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
            % Important Note: the +dG above is correct.
            % The upstream wall moves downstream.
            
            
            %%%need to include choice about circlular perimeter...
            if contains(modelinputs.planshape, 'egg')
                M.xd = M.xd + dC_minor + dE_minor + dM + dG + 0 * dOC +  dP; % if you dont want any melt from dP, then use 0 * dP            
            elseif contains(modelinputs.planshape, 'circle')
                M.xd = M.xd + dC_minor + dE_minor + dM + dG + dOC +  dP; % if you dont want any melt from dP, then use 0 * dP
            end
            
            
            
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
            

            
            
            %%%need to include choice about circlular perimeter...
            if contains(modelinputs.planshape, 'egg')
                M.xu = M.xu - dC_major - 0*dE_major - dM + dG - dOC - 0*dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
            % Important Note: the +dG above is correct.
            % The upstream wall moves downstream.
            
            % The downstream wall also moves downstream
            % at the same rate, dG.
                
                 M.xd = M.xd + dC_minor + 0*dE_minor + dM + dG + 0*dOC +  0*dP; % if you dont want any melt from dP, then use 0 * dP
           
            
            elseif contains(modelinputs.planshape, 'circle')
               M.xu = M.xu - dC_major - 0*dE_major - dM + dG - 0.5*dOC - 0*dP; %melt rate at the apex of the ellipse is 1/2 the total meltrate, which will be nonuniformly distributed along the new perimeter
            % Important Note: the +dG above is correct.
            % The upstream wall moves downstream.
            
            % The downstream wall also moves downstream
            % at the same rate, dG.
                M.xd = M.xd + dC_minor +0* dE_minor + dM + dG + 0.5*dOC +  0*dP; % if you dont want any melt from dP, then use 0 * dP
            end
            
            % Shift them both back upstream so that the bed of the upstream wall
            % stays pinned at x = 0:
            %x0 = M.xu(1);
            %M.xu = M.xu - x0;
            %M.xd = M.xd - x0;
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

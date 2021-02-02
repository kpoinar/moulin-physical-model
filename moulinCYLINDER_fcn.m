% This function does NO geometry evolution.  You start with a cylinder, you
% end with a cylinder.
%
function time = moulinCYLINDER_fcn( workingdirectory, savelocation, makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savefile) 

% 


    cd(workingdirectory)
    
    % Do you want to plot model results?
    make_simple_plots =  makeplots_tf;
    save_figures      =  savefigures_tf;
    visible_figures   =  showfigures_tf;
    
    %% define some basic parameters
    C         = makeConstants;  %constants used for parameterizations  %if there is a particular value we want to change gradually, 
    if isfield(modelinputs,'A_value')%, 'var')
        C.A = modelinputs.A_value; 
        C.c2 = 1*C.A*C.n^(-C.n); % Need to recalculate the closure parameter (Schoof 2010) with the new A
    end
    
    
    Tdatatype = modelinputs.Tdatatype{1}; %'Ryser_foxx';   %ice temperature profile to extrapolate from
    numofdays = modelinputs.numofdays; % 3  %set the number of days for the model run
    H         = modelinputs.H;         % 800   % ice thickness, meters
    R0        = modelinputs.R0;        % 2 m   % radius @bed of moulin initially
    Rtop      = R0;      % 0.5 m % radius @sfc of moulin initially
    L         = modelinputs.L;         % L  % Length of the subglacial channel
    f         = modelinputs.fract_pd_melting;         %f   % fraction of the potential energy used to open the top of the moulin (above water level)
    alpha     = modelinputs.alpha;     % 0.03      % regional surface slope (unitless), for use in Glen's Flow Law
    n         = 0;              % flow law exponent (Glen's Flow Law)
    
    %inital guesses for subglacial model
    hw(1) = H;                  % moulin water level (m)
    S(1)  = 1.5*R0;%1.5* R0;                 % subglacial channel cross sectional area (m^2)
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

        Qbase  = 3.* interp1(baseflow(:,1), baseflow(:,2), time.t, 'spline', 'extrap'); 
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
    

    
    %% save general parameters in time file
    time.parameters.H = H;
    time.parameters.L = L;
    time.parameters.R0 = R0;
    time.parameters.numofdays =  numofdays;
    
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
        Ms   = (pi .* Mrminor_prev .*Mrmajor_prev); %moulin cross-section area
        Mp   = eggperimeter(Mrminor_prev, Mrmajor_prev);   %pi.* (3 .*(Mrminor_prev + Mrmajor_prev) - sqrt((3.* Mrminor_prev + Mrmajor_prev) .* (Mrminor_prev +3 .* Mrmajor_prev))); % wetted/melting perimeter =  ellipse perimeter approx pi [ 3(Mrminor+Mrmajor) - sqrt((3*Mrminor+Mrmajor)(Mrminor+3*Mrmajor))]
        Dh   = (4.*(pi .* Mrminor_prev .* Mrmajor_prev)) ./ Mp; %hydrualic diameter
        Rh   = (pi.* Mrminor_prev .* Mrmajor_prev) ./ Mp; % hydraulic radius
        
        eggp =  eggperimeter(Mrminor_prev, Mrmajor_prev);
        
        %%%%%%%%%%
        % Subglacial Schoof model: Conduit size

        tspan = [0,dt];
        y0    = [hw, S];
        %[hw,S,Qout]   = subglacialsc(Mrminor_prev,z,Qin(cc),H,L,C,tspan,y0);
        opt   = odeset('RelTol', 10.0^(-3), 'AbsTol' , 10.0^(-3));
        %Qin_tot       = Qin(cc) + time.V
        %Qin_compensated = Qin(cc)+Qadd + Qbase(cc); %including Qbase provides a minimum base flow for the subglacial model without needing to route additional water through the moulin
        Qin_compensated = Qin(cc) + Qbase(cc); %including Qbase provides a minimum base flow for the subglacial model without needing to route additional water through the moulin
        time.Qin_compensated(cc) = Qin_compensated;
        [hw,S,Qout, dydt_out]   = subglacialsc(Ms,z,Qin_compensated,H,L,C,dt,tspan,y0, opt); %consider adding Vadd to the qin values
        
            time.S(cc)    = S;
            time.hw(cc)   = hw;
            time.Qout(cc) = Qout ;% - Qbase(cc); %this removes the baseflow from the Qout
            time.dydt_out(cc) = dydt_out;
            time.Qbase(cc) = Qbase(cc); 

    
        %%%%%%%%%%%
        % which nodes are underwater or at the water line (wet) versus above the water line?        
        wet = locatewater(hw,z);
        time.wet(:,cc)    = wet;
        
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
    
    time.dC_minor = zeros(size(time.Mcapacity));
    time.dC_major = zeros(size(time.Mcapacity));
    time.dM = zeros(size(time.Mcapacity));
    time.dP = zeros(size(time.Mcapacity));
    time.dE_minor = zeros(size(time.Mcapacity));
    time.dE_major = zeros(size(time.Mcapacity));
    time.dOC = zeros(size(time.Mcapacity));
    
    %% figures
    
     if make_simple_plots
         simpleplots(time, save_figures, savelocation, visible_figures, modelinputs.runnumber, savefile)
     end
    
    %%

end
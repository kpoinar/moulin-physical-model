function [Qin, Qbase] = Qincalc(Qinreal, Qinfile, Qin_smoothval, Qin_year, Qin_basin, Qin_baseflow, Qin_dampen, ...
                          Qintype, Qinrange, Qinbase, timet, baseflowmultiplier)


% Construct an approximate Qin for each timestep, based on air temps:
    % Qin = double(Tair>C.To) .* (Tair-C.To)*const;
    % Qin = zeros(size(time.t));
    % Use predetermined Qins of various types
    qinreal = Qinreal;

    if qinreal 
        load(Qinfile{1}) %1 = time, 2 cosine function, 3
        Q(:,1) = Q2.(Qin_year{1}).time_seconds;
        Q(:,2) = Qin_dampen .* Q2.(Qin_year{1}).(Qin_basin{1});
        Qin     = interp1(Q(:,1), Q(:,2), timet, 'spline', 'extrap'); % run an interp just in case the timeframe changes
        
        %apply a smoothing to Qin to dampen diurnal varibility
        Qin = max(0.1, Qin);
        Qin = smoothdata(Qin, 'movmean',  Qin_smoothval);
        clear Q
        
        % apply a base flow only to the subglacial system - this kind of
        % mimics a number of processes
        load(Qinfile{1});
        baseflow  =  Q2.(Qin_year{1}).(Qin_baseflow{1});

    
        Qbase  = baseflowmultiplier.* interp1(baseflow(:,1), baseflow(:,2), timet, 'spline', 'extrap'); 

        
        Qbase = Qbase';
        clear Qbase_subglacial baseflow

%lca commented these out in order to prevent Qbase from contributing to moulin geometry          
%         Qin = Qin + Qbase;
%         time.Qin = Qin + Qbase;
        
    else 
        load(Qinfile{1}) %1 = time, 2 cosine function, 3
        Qin     = interp1(Q(:,1), Q(:,Qintype), timet, 'spline', 'extrap'); % run an interp just in case the timeframe changes
        %not sure if we should really have this... might require an adjustment
        %in the qin values
        %Qin     = Qin*modelinputs.Qinmultiplier1 +modelinputs.Qinmultiplier2; %scale Qin to deal with a few model issues
        % Normalize Qin
        Qin = (Qin - mean(Qin)) / range(Qin);
        % Rescale using user-defined Qin parameters
        Qin = max(0.1, Qin * Qinrange + Qinbase);
        time.Qin = Qin;  %save for future plotting
        clear Q 
        
        Qbase = timet .* 0 + 0; %apply zero baseflow to the subglacial model
        
    end
    
figure
hold on 

plot(timet, Qin)
%yyaxis right 
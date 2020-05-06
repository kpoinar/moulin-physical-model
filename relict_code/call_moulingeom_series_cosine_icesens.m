% use this file to run a series of moulin model runs

%% a little set up
clear variables 
close all
clc

savefile = datestr(now, 'mm-dd-yyyy');

if ~exist(['./modeloutputs/', savefile], 'dir')
    disp('making todays file ...')
    mkdir(['./modeloutputs/', savefile]) 
end    

warning('off')
%% variables to define

%set for the entire series of runs
clearpreviousfiles  = true;                 %This clears the files previously saved on this day for this script, use with caution
seriesname          = 'variablethickness';  % a descritive name of the entire runseries
makeplots_tf        = true;              % do you want to make the plots for each?
savefigures_tf      = true;              % do you want to save plots for each?
showfigures_tf      = false;              %show the figures
savelocation        =['./modeloutputs/', savefile];      % where do you want to save the outputs?
workingdirectory    = pwd; %yeah, I know this is a pain, but it makes things work more easily


%set for each rund of the series: use vector for changing values, if a
%single value, then the code at the end will fill in the vector for you
%%%NOTE: variables either need to be 1 or the full number of the runs - I
%%%am not a good enough coder to get 

parameters.H                   = [750];                                     %ice thicknesses to use for the model
                                 [L,alpha, Qinval, Qinmultiplier2 ]  ...
                                     = makeicesheetgeom(parameters.H);      %this function provides a consistent distance from terminus and local slope for
                                                                            % for a given ice thickness based on an idealized ice sheet profile.
%%%For basic sensitivity study, I suggest H = 750, Qin = column 4;
%%%Qinmultiplier1 = 0.75 and Qinmultiplier2 = 4

parameters.numofdays           = 3; %For the sensititivity study, maybe 20d? %For how many days do you want to run the model?

parameters.R0                  = 2;                                          % initial moulin radius



parameters.unfilled_melting    = 1;                                         %what type of melting do you want above the water line? 
                                                                            %1. Open channel, 2. waterfall, 3. potential drop, 4, none
                                                
parameters.A_value             = [6e-24, 6e-24 6e-24 ...
                                  6e-24, 6e-24, 6e-24 6e-24 6e-24 ...
                                  1e-24 5e-24 10e-24]   ;                                %A value for the subglacial
                                                                            %channel - only uncomment if you want to change this with each model run,
                                                                            %otherwise, set in makeConstants
parameters.f_sub               = 0.05;                                      % friction factor for the subglacial system - it should be ~ between 0.03 and 0.2, the water level in the moulin is very sensitive to this.
                                                
parameters.Tdatatype           = {'Ryser_foxx','Luthi', 'HarrS4C' ... 
                                 'Ryser_foxx','Ryser_foxx','Ryser_foxx','Ryser_foxx' ...
                                 'Ryser_foxx','Ryser_foxx','Ryser_foxx'};
                                                                            % Choose the ice temperature to use for the model  
                                
parameters.L                   = [L];                                       % distance from terminus for subglaical channel

parameters.f                   = 0.1;                                       % this is the fraction of melting for potential drop 

parameters.alpha               = [alpha];                                   % regional slope

parameters.E                   = [5 5 5 ...
                                  1 3 5 7 10 ...
                                  5 5 5];                                % moulin deformation enhancement factor for creep

parameters.Qinfile             = {'Q_cosine_elevation_bands.mat'};          % to use Qin value 
%parameters.Qinfile             = {'Qcosines_new.mat'};                     %This is the older file 
parameters.Qintype             = [4];%Qinval                                % Qin type: 1. cosine, 2. cosine with small melt events, 
                                                                            % 3. cosine with large melt events, 4. quasi-real data, 
                                                                            % 5. realistic, 6 realistic with taper to near zero
                                                   
parameters.Qinmultiplier1      = 0.75; %1                                   % of the form: Qin     = Qin*Qinmultiplier1 + Qinmultiplier2;
parameters.Qinmultiplier2      = 4 ; %Qinmultiplier2                        % of the form: Qin     = Qin*Qinmultiplier1 + Qinmultiplier2; %This needs to be tuned a litle bit

parameters.relative_roughness  = 0.2;                                       %relative roughness for below the water line, increasing will increase the amount of melt due to turbulence

parameters.relative_roughness_OC =  0.001;  %0.1;                                    %relative roughness above the water line for open channel

parameters.fR_wet                = [1];                                     % This is used when in moulingeom_fcn, fR_fixed is true (right now it is on 4/7/20 LCA)
%% find the longest series and repeat the other ones and do some basic file maintainence

savetime = datestr(now, 'mm-dd-yyyy_HHMM');
maxlengthofseries    = max(structfun(@length,parameters));
parameters.runnumber = 1:maxlengthofseries;
lengthofseries       = structfun(@length,parameters);
names                = fieldnames(parameters);

for ii = 1:length(lengthofseries)
    if lengthofseries(ii) < maxlengthofseries
        parameters.(names{ii}) = repmat(parameters.(names{ii}), 1, maxlengthofseries);
    end
end


if clearpreviousfiles %to be safe rm is in interactive mode
    cd(savelocation)
    !rm -i *   
end
    
%%

for ii = 1:maxlengthofseries
    cd(workingdirectory)
    disp(['Working on run ', num2str(ii), ' of ' num2str(maxlengthofseries), '...'])
    for jj = 1:length(lengthofseries)
        modelinputs.(names{jj}) = parameters.(names{jj})(ii);
    end
    time   =  moulingeom_fcn(workingdirectory, savelocation,  makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savetime);
    seriesresults{ii,1} = ['run_', num2str(ii)];
    seriesresults{ii,2} = time;
    seriesresults{ii,3} = modelinputs;
    %clear time
end




   disp(['Saving ', seriesname, ' in ', savelocation,  '...'])
   cd(savelocation)

   filename = ['modelrunseries','-',seriesname,  '-', savetime, '_outputs.mat'];
   save(filename, 'seriesresults')
   clearvars -except seriesresults workingdirectory savelocation savefile
   
   %copy script with new name and tar entire package
   copyfile([workingdirectory, '/call_moulingeom_series_cosine_icesens.m'], [savelocation, '/', savefile, '_call_moulingeom_series_cosine_ice_sens.m']);
   
   disp('Done!')


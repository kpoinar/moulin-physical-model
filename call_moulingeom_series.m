% use this file to run a series of moulin model runs

clear variables 
close all
clc

savefile = datestr(now, 'mm-dd-yyyy');

if ~exist(['~/Documents/Repositories/Moulin/moulin-physical-model/modeloutputs/', savefile], 'dir')
    disp('making todays file ...')
    mkdir(['~/Documents/Repositories/Moulin/moulin-physical-model/modeloutputs/', savefile]) 
end    

warning('off')
%% variables to define

%set for the entire series of runs
clearpreviousfiles  = true;                 %This clears the files previously saved on this day for this script, use with caution
seriesname          = 'variablethickness';  % a descritive name of the entire runseries
makeplots_tf        = true;              % do you want to make the plots for each?
savefigures_tf      = true;              % do you want to save plots for each?
showfigures_tf      = false;              %show the figures
savelocation        =['~/Documents/Repositories/Moulin/moulin-physical-model/modeloutputs/', savefile];      % where do you want to save the outputs?
workingdirectory    = '~/Documents/Repositories/Moulin/moulin-physical-model'; %yeah, I know this is a pain, but it makes things work more easily



%set for each rund of the series: use vector for changing values, if a
%single value, then the code at the end will fill in the vector for you
%%%NOTE: variables either need to be 1 or the full number of the runs - I
%%%am not a good enough coder to get 

parameters.unfilled_melting    = 1;              %what type of melting do you want above the water line? 
                                                  %1. Open channel, 2. waterfall, 3. potential drop, 4, none
                                                
%parameters.A_value             = nan;           %A value for the subglacial
                                                  %channel - only uncomment if you want to change this with each model run,
                                                  %otherwise, set in makeConstants
                                                
parameters.Tdatatype           = {'Ryser_foxx'};   %Choose the ice temperature to use for the model                                                

parameters.numofdays           = 5;              %for how many days do you want to run the model?

parameters.H                   = [800 500 1200];  %ice thicknesses to use for the model

parameters.R0                  = 1;               % initial moulin radius

parameters.L                   = [12e3 8e3 15e3]; % distance from terminus for subglaical channel

parameters.f                   = 0.07;            %subglacial friction factor

parameters.alpha               = 0.03;            %regional slope

parameters.E                   = 5;               %moulin deformation enhancement factor for creep

parameters.Qintype             = 5;               %Qin type: 1. cosine, 2. cosine with small melt events, 
                                                   % 3. cosine with large melt events, 4. quasi-real data, 
                                                   % 5. realistic, 6 realistic with taper to near zero
                                                   
parameters.Qinmultiplier1      = 0.8;             % of the form: Qin     = Qin*Qinmultiplier1 + Qinmultiplier2;
parameters.Qinmultiplier2      = 3;               % of the form: Qin     = Qin*Qinmultiplier1 + Qinmultiplier2;

parameters.relative_roughness  = 0.2;             %relative roughness for below the water line, increasing will increase the amount of melt due to turbulence

parameters.relative_roughness_OC = 1e-9;          %relative roughness above the water line for open channel

%% find the longest series and repeat the other ones...

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
   clearvars -except seriesresults
   disp('Done!')


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


parameters.numofdays           = 60;              %for how many days do you want to run the model?


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

parameters.H                   = fliplr([669 820 947 1058 1159 1252 1339 1420 1497 1569]);  %ice thicknesses to use for the model


[L,alpha, Qinrange, Qinbase ]  = makeicesheetgeom(parameters.H);      %this function provides a consistent distance from terminus and local slope for
                                                   % for a given ice thickness based on an idealized ice sheet profile.

parameters.R0                  = 2;               % initial moulin radius

parameters.L                   = [L];%[20e3 30e3 15e3]; % distance from terminus for subglaical channel

parameters.f                   = 0.07;            %subglacial friction factor

parameters.alpha               = [alpha];% 0.03;            %regional slope

parameters.E                   = 5;               %moulin deformation enhancement factor for creep

parameters.Qinfile             = {'Q_cosine_elevation_bands.mat'};

parameters.Qintype             = 2;               % 
                                                   
parameters.Qinrange            = 0.5;     % of the form:  Qin     = Qin*Qinrange + Qinbase;
%     nqinrange = length(parameters.Qinrange);
parameters.Qinbase             = 3.5;      % of the form:  Qin     = Qin*Qinrange + Qinbase;
%     nqinbase  = length(parameters.Qinbase);

% temp = parameters.Qinrange;
% parameters.Qinrange  = repmat(parameters.Qinrange,1,size(parameters.Qinbase));
% parameters.Qinbase   = sort(repmat(parameters.Qinbase,length(temp),1)'); clear temp

parameters.relative_roughness  = 0.2;             %relative roughness for below the water line, increasing will increase the amount of melt due to turbulence

parameters.relative_roughness_OC = 1e-9;          %relative roughness above the water line for open channel

parameters.fR_wet                = 0.1   
%% find the longest series and repeat the other ones and do some basic file maintainence

savetime = datestr(now, 'mm-dd-yyyy_HHMM');
maxlengthofseries    = max(structfun(@length,parameters));
parameters.runnumber = 1:maxlengthofseries;
lengthofseries       = structfun(@length,parameters);
names                = fieldnames(parameters);

tic
for ii = 1:length(lengthofseries)
    if lengthofseries(ii) < maxlengthofseries
        parameters.(names{ii}) = repmat(parameters.(names{ii}), 1, maxlengthofseries);
    end
end


if clearpreviousfiles %to be safe rm is in interactive mode
    cd(savelocation)
    %!rm -i *   
    !rm *.png
    cd(workingdirectory)
end
    
%%

for ii = 1:maxlengthofseries
    cd(workingdirectory)
%     disp(['Working on run ', num2str(ii), ' of ' num2str(maxlengthofseries), '...'])
    fprintf('\nWorking on run %d of %d:   Qin = %1.1f range + %1.1f base   at H=%d m: \n',ii,maxlengthofseries,parameters.Qinrange(ii),parameters.Qinbase(ii),parameters.H(ii));
    for jj = 1:length(lengthofseries)
        modelinputs.(names{jj}) = parameters.(names{jj})(ii);
    end
    time   =  moulingeom_fcn(workingdirectory, savelocation,  makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savetime);
    seriesresults{ii,1} = ['run_', num2str(ii)];
    seriesresults{ii,2} = time;
    seriesresults{ii,3} = modelinputs;
    
end

%
%
%%

   disp(['Saving ', seriesname, ' in ', savelocation,  '...'])
   cd(workingdirectory)
   cd(savelocation)

   filename = ['modelrunseries','-',seriesname,  '-', savetime, '_outputs.mat'];
   save(filename, 'seriesresults')   
   clearvars -except seriesresults workingdirectory savelocation savefile success
   cd(workingdirectory)
   
   %copy script with new name and tar entire package
   copyfile([workingdirectory, '/call_moulingeom_series.m'], [savelocation, '/', savefile, '_call_moulingeom_series_cosine_watersens.m']);
   
   disp('Done!')
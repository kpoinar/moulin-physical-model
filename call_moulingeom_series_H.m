% use this file to run a series of moulin model runs

%% a little set up
clear variables 
%close all
%clc

seriesname          = 'real';  % a descritive name of the entire runseries
savefile = datestr(now, 'mm-dd-yyyy');

if ~exist(['./modeloutputs/', seriesname, savefile], 'dir')
    disp('making todays file ...')
    mkdir(['./modeloutputs/', seriesname, savefile]) 
    %mkdir(['~/Documents/Repositories/Moulin/real/moulin-physical-model/modeloutputs/', savefile]) 
end    

warning('off')
%% variables to define

%set for the entire series of runs
clearpreviousfiles  = false;                 %This clears the files previously saved on this day for this script, use with caution
makeplots_tf        = true;              % do you want to make the plots for each?
savefigures_tf      = false;              % do you want to save plots for each?
showfigures_tf      = true;              %show the figures
savelocation        =['./modeloutputs/' ,seriesname, savefile];
% savelocation         = ['~/Documents/Repositories/Moulin/real/moulin-physical-model/modeloutputs/', savefile]; % where do you want to save the outputs?
workingdirectory    = pwd; %yeah, I know this is a pain, but it makes things work more easily


parameters.numofdays           =  10 ;              %for how many days do you want to run the model?


%set for each rund of the series: use vector for changing values, if a
%single value, then the code at the end will fill in the vector for you
%%%NOTE: variables either need to be 1 or the full number of the runs - I
%%%am not a good enough coder to get 

parameters.unfilled_melting    = 1;              %what type of melting do you want above the water line? 
                                                  %1. Open channel, 2. waterfall, 3. potential drop, 4, none
                                                
parameters.A_value             = [ 6e-24];           %A value for the subglacial
                                                  %channel - only uncomment if you want to change this with each model run,
                                                  %otherwise, set in makeConstants
                                                  
%parameters.ShearModulus        = [1  3  5  7  9]*1e9;     % Shear modulus (elastic behavior)
                                                  %only uncomment if you want to change this with each model run,
                                                  %otherwise, set in makeConstants
                                                
parameters.Tdatatype           = {'Ryser_foxx'};   %Choose the ice temperature to use for the model                                                
%parameters.Tdatatype  = {'IkenB','Ryser_gull','Ryser_foxx','Cold','Temperate'};

% parameters.H                   = ([669 820 947 1058 1159 1252 1339 1420 1497 1569]);  % 10 ice thicknesses - may work on CCR
%parameters.H                   = ([669  820  1058  1339  1569]);  % 5 ice thicknesses (kulusuk limit for 10 day runs)
parameters.H                   = 820;    % H = 820 is our standard moulin, 30 km from the margin.


% ID	surf	bed	thickness
% 44	949	    396	553**** done
% 184	1254	240	1014
% 352	1144	403	741 *** done 
% 718	1623	203	1420
% 722	1554	239	1315**
% 727	1307	130	1177
% 733	1363	278	1085
% 749	978	    -35	1013  **
% 754	1051	154	897
% 895	1514	225	1289
% 911	1539	200	1339


[L,alpha, Qinrange, Qinbase ]  = makeicesheetgeom(parameters.H);      %this function provides a consistent distance from terminus and local slope for
                                                   % for a given ice thickness based on an idealized ice sheet profile.
                                                   
%parameters.R0            = [0.65  1.0   2.0  3.5  5.0];
parameters.R0                  = 2;%0.75;               % initial moulin radius @ bed
parameters.Rtop                = 0.5;             % initial moulin radius @ surface (should be less than R0)

parameters.L                   = L;%[20e3 30e3 15e3]; % distance from terminus for subglaical channel

parameters.use_pD_downstream   = false;%true; 
parameters.fract_pd_melting    = 0.0;            %This is fraction of energy from potential drop that is used for melting
%This is currently off in the model



parameters.alpha               = [alpha];% 0.03;            %regional slope

parameters.E                   = [5];               %moulin deformation enhancement factor for creep



%%%%%%%%%%%
parameters.Qinreal             = false;
% if Qin is from realistic data
if parameters.Qinreal
    parameters.Qinfile             = {'/Users/lcandre2/Documents/Repositories/Moulin/real/moulin-physical-model/Qin_real_2019.mat'}; %{'Qincosines_peak7pm.mat'}; '/
    
    parameters.Qin_smoothval       = 12; % for the smooth function assuming 1 hour data. This dampens diurnal varibility
    parameters.Qin_year            = {'y2019'};
    parameters.Qin_basin           = {'b722'}; %'b749' %b727
    parameters.Qin_baseflow        = {'b722_baseflow'};  %basename_baseflow, column 1 time, column2 baseflow estimate
    parameters.Qin_dampen          = 0.5; %this dampens the diurnal varibility
    %%%%%%%%%%
else
    
    % if Qin is idealized
    parameters.Qinfile             = {'Q_cosine_elevation_bands.mat'} %{'Qincosines_peak7pm.mat'}; '/
    parameters.Qintype             = 2;    % 2 is cosinusoid, 3 is two superimposed cosinusoids, 4 has more variability, 5 has most variability
    parameters.Qinrange            = 0.5;     % of the form:  Qin     = Qin*Qinrange + Qinbase;
    %     nqinrange = length(parameters.Qinrange);
    parameters.Qinbase             = 3.5;      % of the form:  Qin     = Qin*Qinrange + Qinbase;
    %     nqinbase  = length(parameters.Qinbase);
end
%%%%%%%%%%

% temp = parameters.Qinrange;
% parameters.Qinrange  = repmat(parameters.Qinrange,1,size(parameters.Qinbase));
% parameters.Qinbase   = sort(repmat(parameters.Qinbase,length(temp),1)'); clear temp

parameters.variable_fR_oc        = false;
parameters.variable_fR_wet       = false;
%if variable_fR_oc and variable_fR_wet are false, then these do *NOT* matter
parameters.relative_roughness    = 0.2;             %relative roughness for below the water line, increasing will increase the amount of melt due to turbulence
parameters.relative_roughness_OC = 0.5;          %relative roughness above the water line for open channel
 
%if variable_fR_oc and variable_fR_wet are true, then these do *NOT* matter
parameters.fR_fixed_wet                = 0.1;%logspace(-2,0,5);%0.1;  
parameters.fR_fixed_oc                 = 0.5; %logspace(-2,0,5); %1;  % base value was 10!

% Background stress conditions (for elastic deformation)
parameters.sigx = 0e3;    % compressive (-) or extensive (+).  Make sigy + and keep sigx 0.
parameters.sigy = 50e3;    % compressive (-) or extensive (+).  Make sigy + and keep sigx 0.  
parameters.tauxy = -50e3;  % shear opening; compressive (-) or extensive (+)
% sigx-sigy is about -7x more important than tauxy in determining the surface
% expression of the moulin
% we want net stress at sfc to be compressive, so sigx-sigy (+) and/or tauxy (-)
% Strain rate data (Figure 1 of "Challenges" paper) indicate about -30 kPa
% mean principal stress at moulin sites


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
    !rm *   
    %y!rm *.mat
    cd(workingdirectory)
end
    
%%
% cell for results
res = cell(maxlengthofseries,3);

for ii = 1:maxlengthofseries
    cd(workingdirectory)
     disp(['Working on run ', num2str(ii), ' of ' num2str(maxlengthofseries), '...'])
    %fprintf('\nWorking on run %d of %d:   at H=%d m: \n',ii,maxlengthofseries,parameters.Qinrange(ii),parameters.Qinbase(ii),parameters.H(ii));
    for jj = 1:length(lengthofseries)
        modelinputs.(names{jj}) = parameters.(names{jj})(ii);
    end
    time   =  moulingeom_fcn(workingdirectory, savelocation,  makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savetime);
   
    % DO STATIC CYLINDER MOULIN
    %time   =  moulinCYLINDER_fcn(workingdirectory, savelocation,  makeplots_tf, savefigures_tf, showfigures_tf, modelinputs, savetime);

    
    res{ii,1} = ['run_', num2str(ii)];
    res{ii,2} = time;
    res{ii,3} = modelinputs;
    
end

%
%
%%

   disp(['Saving ', seriesname, ' in ', savelocation,  '...'])
   cd(workingdirectory)
   %cd(savelocation)

   filename = [savelocation, '/modelrunseries','-',seriesname,  '-', savetime, '_outputs.mat'];
   save(filename, 'res')   
   clearvars -except res workingdirectory savelocation savefile success
   cd(workingdirectory)
   
   %copy script with new name and tar entire package
   copyfile([workingdirectory, '/call_moulingeom_series_H.m'], [savelocation, '/', savefile, '_call_moulingeom_series_used.m']);
   
   disp('Done!')
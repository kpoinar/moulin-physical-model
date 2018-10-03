clear all
close all
clc

cd ~/Documents/Projects/Quiescent/reduce/Moulin_model/Plotting/Melt/
load foxx_temp_melt_rate_2011.mat
load gull_temp_melt_rate_2011.mat
load hare_temp_melt_rate_2011.mat

%% Load the merra-2 data and create a timeseries of the x,y, runoff for inputs

cd ~/Documents/Projects/moulin_formation/data/m2_runoff/
files        = dir('*2016*.nc4');
files        = struct2cell(files);
files        = files(1, 1:end);
files   = files';

%%
m2_runoff = nan(109,60);
m2_time   = nan;
for ii = 1:length(files)
    tmp(:,:,:) = ncread(([files{ii}]), 'RUNOFF');
    m2_runoff  = cat(3, m2_runoff, tmp);
    
    tmp2       = ncread(([files{ii}]), 'time');
    m2_time    = cat(1, m2_time, tmp2);
    year(ii,1)   = datenum(str2num(files{ii}(28:31)));
    month(ii,1)   = datenum(str2num(files{ii}(32:33)));
    day(ii,1)     = datenum(str2num(files{ii}(34:35)));
end

m2_runoff = m2_runoff(:,:,2:end);

%%
lat = ncread(([files{1}]), 'lat');
lon = ncread(([files{1}]), 'lon');

%[lat, lon] = meshgrat(lat, lon);
% lat = lat'; lon = lon';
% %%
% 
dateall = datenum(year,month, day);
dateall = repmat(dateall, 8,1);
dateall = sort(dateall);
dateall = dateall + (m2_time(2:end) +90)/1440;
% 
% 
% %% find coordinates of boreholes for foxx and russell
% %zwally_outline = importdata('~/Documents/Data/greenland_basinz/zwally_outline.txt');
% %zlat = zwally_outline.data(:,2);
% %zlon = zwally_outline.data(:,3);
% %zwally_basins  = importdata('~/Documents/Data/greenland_basinz/zwally_basins.txt');
% %zbasin_lat = zwally_basins.data(:,2);
% %zbasin_lon = zwally_basins.data(:,3);
% 
% %load quantarctica basemap outline 
% %quant = shaperead('~/Documents/Data/Quantarctica2/Basemap/World/NE_50m_Countries.shp');
% 
% fractice = ncread('~/Documents/Data/MERRA-2/constant_parameters/MERRA2_101.const_2d_asm_Nx.00000000.nc4','FRLANDICE' );
% fractice = fractice(169:277, 298:357) ;
% 
% %% 
% 
% harrington = m2_runoff(41,19,:); 
% %%
% 
% harr_fract = fractice(41,18)
% 
% 
% 
% %%
% load coastlines
% 
%  figure
%  hold on; box on; grid on
%   axesm ('lambertstd','MapLatLimit',[59 84],'MapLonLimit',[-75 0])
%   axis off; framem on; gridm on; mlabel on; plabel on;
%    geoshow(coastlat,coastlon,'DisplayType','line', 'color', 'k', 'linewidth', 1);
% 
%  plotm(lat(41,19), lon(41,19), '*r');



%%
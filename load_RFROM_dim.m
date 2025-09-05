% load_RFROM_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of the RFROM 
% temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/23/2024

function [TS,months,weeks,timesteps] = load_RFROM_dim(fpath,y1,y2)

% % check for existence of file
% if ~isfile([fpath 'RFROMV22_TEMP_STABLE_CLIM_' num2str(y1) '_' num2str(y2) '.nc'])
% 
% end

% load RFROM climatological temp and salinity
TS.Longitude = ncread([fpath 'RFROMV22_TEMP_STABLE_CLIM_' num2str(y1) '_' num2str(y2) '.nc'],'longitude');
TS.Latitude = ncread([fpath 'RFROMV22_TEMP_STABLE_CLIM_' num2str(y1) '_' num2str(y2) '.nc'],'latitude');
TS.Pressure = ncread([fpath 'RFROMV22_TEMP_STABLE_CLIM_' num2str(y1) '_' num2str(y2) '.nc'],'mean_pressure');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Pressure);
% determine number of monthly timesteps
files = dir([fpath '/RFROM_TEMP_v2.2/*.nc']);
for t = 1:length(files)
    file_date = char(extractBetween(files(t).name,'STABLE_','.nc'));
    file_year(t) = str2num(file_date(1:4));
end
idx1 = find(file_year == y1,1,'first');
idx2 = find(file_year == y2,1,'last');
% process time
months = [];
weeks = [];
timesteps = 0;
TS.Time = [];
TS.cnt = {};
cnt = 1;
for m = idx1:idx2
    % determine number of weeks in each monthly file
    nc_atts = ncinfo([files(m).folder '/' files(m).name]);
    num_weeks = nc_atts.Dimensions(3).Length;
    months = [months;repmat((m+1)-idx1,num_weeks,1)];
    weeks = [weeks;(1:num_weeks)'];
    timesteps = timesteps + num_weeks;
    TS.Time = [TS.Time;double(ncread([files(m).folder '/' files(m).name],'time'))];
    TS.cnt = [TS.cnt;cnt:(cnt-1)+num_weeks];
    cnt = cnt + num_weeks;
end
% add years and months
TS.years = nan(length(files),1);
TS.months = nan(length(files),1);
for n = 1:length(files)
    date = cell2mat(extractBetween(files(n).name,'STABLE_','.nc'));
    TS.years(n) = str2double(date(1:4));
    TS.months(n) = str2double(date(6:7));
end
% remove months and years outside range
idx = TS.years < y1 | TS.years > y2;
TS.years(idx) = [];
TS.months(idx) = [];

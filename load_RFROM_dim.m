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

function [TS,timesteps] = load_RFROM_dim(fpath)

% load RFROM climatological temp and salinity
TS.Longitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'longitude');
TS.Latitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'latitude');
TS.Pressure = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'mean_pressure');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Pressure);
% determine number of monthly timesteps
files = dir([fpath 'RFROM_TEMP_v0.1/*.nc']);
timesteps = length(files);
% process time
TS.years = nan(length(files),1);
TS.months = nan(length(files),1);
for n = 1:length(files)
    date = cell2mat(extractBetween(files(n).name,'STABLE_','.nc'));
    TS.years(n) = str2double(date(1:4));
    TS.months(n) = str2double(date(6:7));
end
idx = TS.years < 2004;
TS.years(idx) = [];
TS.months(idx) = [];

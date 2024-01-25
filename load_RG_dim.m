% load_RG_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of the Roemmich 
% and Gilson temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/23/2024

function [TS,timesteps] = load_RG_dim(fpath)

% load temperature and salinity
TS.Longitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Longitude');
TS.Latitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Latitude');
TS.Pressure = ncread([fpath 'RG_Climatology_Temp.nc'],'Pressure');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Pressure);
% determine number of monthly timesteps
TS.Time = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time');
timesteps = length(TS.Time);
% process time
dates = datevec(datenum(2004,1,1+double(TS.Time)));
TS.years = dates(:,1);
TS.months = dates(:,2);
clear dates

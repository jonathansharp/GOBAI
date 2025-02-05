% load_model_dim
%
% DESCRIPTION:
% This function is used to load the 
% spatiotemporal dimensions of CMIP
% model output.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/1/2024

function [TS,months,weeks,timesteps] = load_model_dim(fpath)

% load dimensions
TS.Latitude = ncread(fpath,'lat');
TS.Longitude = ncread(fpath,'lon');
TS.Depth = ncread(fpath,'depth');
% compute dimensions
TS.xdim = length(TS.Longitude);
TS.ydim = length(TS.Latitude);
TS.zdim = length(TS.Depth);
% determine number of monthly timesteps
TS.Time = ncread(fpath,'time');
timesteps = length(TS.Time);
months = 1:timesteps;
weeks = ones(timesteps,1);
% add years and months
dates = datevec(datenum(0,0,1+double(TS.Time)));
TS.years = dates(:,1);
TS.months = dates(:,2);
clear dates

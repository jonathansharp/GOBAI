% load_RFROM_clim
%
% DESCRIPTION:
% This function is used to load climatological 
% temperature and salinity from the RFROM 
% temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/23/2024

function TS = load_RFROM_clim(TS,fpath)

% load RFROM climatological temp and salinity
TS.temperature_cns = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'ocean_temperature');
TS.salinity_abs = ncread([fpath 'RFROM_SAL_STABLE_CLIM.nc'],'ocean_salinity');
TS.temperature_cns = mean(TS.temperature_cns,4,'omitnan');
TS.salinity_abs = mean(TS.salinity_abs,4,'omitnan');

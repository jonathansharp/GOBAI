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
TS.temp_clim = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'ocean_temperature');
TS.sal_clim = ncread([fpath 'RFROM_SAL_STABLE_CLIM.nc'],'ocean_salinity');
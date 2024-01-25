% load_RG_clim
%
% DESCRIPTION:
% This function is used to load climatological 
% temperature and salinity from the Roemmich 
% and Gilson temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/23/2024

function TS = load_RG_clim(TS,fpath)

% load temperature and salinity
TS.Temperature = ncread([fpath 'RG_Climatology_Temp.nc'],'Temperature');
TS.Salinity = ncread([fpath 'RG_Climatology_Sal.nc'],'Salinity');
% calculate climatological mean temperature and salinity
TS.temperature = single(nan(TS.xdim,TS.ydim,TS.zdim,12));
TS.salinity = single(nan(TS.xdim,TS.ydim,TS.zdim,12));
for m = 1:12
    TS.temperature(:,:,:,m) = mean(TS.Temperature(:,:,:,m:12:end),4,'omitnan');
    TS.salinity(:,:,:,m) = mean(TS.Salinity(:,:,:,m:12:end),4,'omitnan');
end
% clean up
TS = rmfield(TS,{'Temperature' 'Salinity'});
% load_model_clim
%
% DESCRIPTION:
% This function is used to load climatological 
% temperature and salinity from CMIP model output.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/1/2024

function TS = load_model_clim(TS,nc_filepath_abs_sal,nc_filepath_cns_tmp)

% load temperature and salinity
TS.abs_sal = ncread(nc_filepath_abs_sal,'abs_sal');
TS.cns_tmp = ncread(nc_filepath_cns_tmp,'cns_tmp');
% calculate climatological mean temperature and salinity
TS.salinity_abs = single(nan(TS.xdim,TS.ydim,TS.zdim,12));
TS.temperature_cns = single(nan(TS.xdim,TS.ydim,TS.zdim,12));
for m = 1:12
    TS.salinity_abs(:,:,:,m) = mean(TS.abs_sal(:,:,:,m:12:end),4,'omitnan');
    TS.temperature_cns(:,:,:,m) = mean(TS.cns_tmp(:,:,:,m:12:end),4,'omitnan');
end
% clean up
TS = rmfield(TS,{'cns_tmp' 'abs_sal'});
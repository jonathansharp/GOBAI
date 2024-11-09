% load_climatological_data
% 
% DESCRIPTION:
% This function loads climatological data
% from a given gridded dataset
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/18/2024

function [TS,timesteps] = load_climatological_TS_data(fpath,base_grid,start_year,date_str)

%% load climatological data
if strcmp(base_grid,'RG')
    % load practical salinity and in situ temperature from RG grid
    [TS,timesteps] = load_RG_dim([fpath '/Data/RG_CLIM/']);
    TS = load_RG_clim(TS,[fpath '/Data/RG_CLIM/']);
    % calculate absolute salinity and conservative temperature
    TS = replicate_dims(base_grid,TS,12);
    idx = ~isnan(TS.salinity) & ~isnan(TS.temperature);
    TS.salinity_abs = nan(size(idx));
    TS.salinity_abs(idx) = gsw_SA_from_SP(TS.salinity(idx),TS.pressure(idx),...
        convert_lon(TS.longitude(idx)),TS.latitude(idx));
    TS.temperature_cns = nan(size(idx));
    TS.temperature_cns(idx) = ...
        gsw_CT_from_t(TS.salinity_abs(idx),TS.temperature(idx),TS.pressure(idx));
    TS = rmfield(TS,{'temperature' 'salinity'});
    % load chlorophyll
elseif strcmp(base_grid,'RFROM')
    % load absolute salinity and conservative temperature from RFROM grid
    [TS,timesteps] = load_RFROM_dim([fpath '/Data/RFROM/']);
    TS = load_RFROM_clim(TS,[fpath '/Data/RFROM/']);
    TS = replicate_dims(base_grid,TS,1);
    % load chlorophyll 
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = '_r1i1p1f1_gr';
    % define filepaths
    nc_filepath_abs_sal = [fpath 'combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    nc_filepath_cns_tmp = [fpath 'combined/regridded/cns_tmp' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load absolute salinity and conservative temperature from model output
    [TS,timesteps] = load_model_dim(nc_filepath_abs_sal);
    TS = load_model_clim(TS,nc_filepath_abs_sal,nc_filepath_cns_tmp);
    TS = replicate_dims(base_grid,TS,12);
end

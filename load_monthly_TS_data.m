% load_climatological_data
% 
% DESCRIPTION:
% This function loads climatological data
% from a given gridded dataset
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/18/2024

function TS = load_monthly_TS_data(fpath,base_grid,m,w,start_year,date_str)

%% load monthly data
if strcmp(base_grid,'RG')
    % load dimensions
    TS = load_RG_dim([fpath '/Data/RG_CLIM/']);
    TS = replicate_dims(base_grid,TS,1);
    % get RG practical salinity and in situ temperature
    TS.temperature = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],...
        'Temperature',[1 1 1 m],[Inf Inf Inf 1]);
    TS.salinity = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Sal.nc'],...
        'Salinity',[1 1 1 m],[Inf Inf Inf 1]);
    % convert to absolute salinity and conservative temperature
    idx = ~isnan(TS.salinity) & ~isnan(TS.temperature);
    TS.salinity_abs = nan(size(idx));
    TS.salinity_abs(idx) = gsw_SA_from_SP(TS.salinity(idx),TS.pressure(idx),...
        convert_lon(TS.longitude(idx)),TS.latitude(idx));
    TS.temperature_cns = nan(size(idx));
    TS.temperature_cns(idx) = ...
        gsw_CT_from_t(TS.salinity_abs(idx),TS.temperature(idx),TS.pressure(idx));
elseif strcmp(base_grid,'RFROM')
    % load dimensions
    TS = load_RFROM_dim([fpath '/Data/RFROM/']);
    TS = replicate_dims(base_grid,TS,1);
    % get RFROM absolute salinity and conservative temperature
    TS.temperature_cns = ncread([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
    num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
        'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
    TS.salinity_abs = ncread([fpath '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
    num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
        'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = '_r1i1p1f1_gr';
    % define filepaths
    nc_filepath_abs_sal = [fpath 'combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    nc_filepath_cns_tmp = [fpath 'combined/regridded/cns_tmp' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
    TS = replicate_dims(base_grid,TS,1);
    % get RG practical salinity and in situ temperature
    TS.salinity_abs = ncread(nc_filepath_abs_sal,'abs_sal',[1 1 1 m],[Inf Inf Inf 1]);
    TS.temperature_cns = ncread(nc_filepath_cns_tmp,'cns_tmp',[1 1 1 m],[Inf Inf Inf 1]);
end

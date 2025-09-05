function fpaths = path_config(base_grid,param)

% model path
fpaths.model_path = '/fast7/model/';

% temperature and salinity paths
if strcmp(base_grid,'RG')
    fpaths.temp_path = [pwd '/Data/RG_CLIM/'];
    fpaths.sal_path = [pwd '/Data/RG_CLIM/'];
elseif strcmp(base_grid,'RFROM')
    fpaths.temp_path = '/fast2/temperature/';
    fpaths.sal_path = '/fast3/salinity/';
else
    fpaths.temp_path = '';
    fpaths.sal_path = '';
end

% parameter path
if strcmp(param,'o2')
    fpaths.param_path_temp = '/fast4/o2/';
    fpaths.param_path = '/data2/o2/';
elseif strcmp(param,'no3')
    fpaths.param_path_temp = '/fast5/no3/';
    fpaths.param_path = '/data2/no3/';
elseif strcmp(param,'dic')
    fpaths.param_path_temp = '/fast6/dic/';
    fpaths.param_path = '/data2/dic/';
end

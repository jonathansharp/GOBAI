function [model_path,param_path,temp_path,sal_path] = path_config(base_grid,param)

% model path
model_path = '/fast7/model/';

% temperature and salinity paths
if strcmp(base_grid,'RG')
    temp_path = [pwd '/Data/RG_CLIM/'];
    sal_path = [pwd '/Data/RG_CLIM/'];
elseif strcmp(base_grid,'RFROM')
    temp_path = '/fast2/temperature/';
    sal_path = '/fast3/salinity/';
else
    temp_path = '';
    sal_path = '';
end

% parameter path
if strcmp(param,'o2')
    param_path = '/fast4/o2/';
elseif strcmp(param,'no3')
    param_path = '/fast5/no3/';
elseif strcmp(param,'dic')
    param_path = '/fast6/dic/';
end

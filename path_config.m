function [model_path,param_path,temp_path,sal_path] = path_config(base_grid,param)

% model path
model_path = '/fast7/sharp/model/';

% temperature and salinity paths
if strcmp(base_grid,'RG')
    temp_path = [pwd '/Data/RG_CLIM/'];
    sal_path = [pwd '/Data/RG_CLIM/'];
elseif strcmp(base_grid,'RFROM')
    temp_path = '/fast2/sharp/temp/';
    sal_path = '/fast3/sharp/sal/';
else
    temp_path = '';
    sal_path = '';
end

% parameter path
if strcmp(param,'o2')
    param_path = '/fast4/sharp/o2/';
elseif strcmp(param,'no3')
    param_path = '/fast5/sharp/no3/';
elseif strcmp(param,'dic')
    param_path = '/fast6/sharp/dic/';
end

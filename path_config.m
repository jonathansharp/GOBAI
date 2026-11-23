function [fpaths,param_props] = path_config(base_grid,param,param_props)

if nargin < 3
    param_props = struct();
end

% shared external paths
fpaths.ancillary_data_path = '/med2/gobai_ancillary_data/';
fpaths.snapshot_dir = '/med2/';
fpaths.snapshot_path = '/med2/BGC_Argo_Snapshots/';
fpaths.figures_path = '/data2/GOBAI/Figures/';
param_props.ancillary_data_path = fpaths.ancillary_data_path;
param_props.snapshot_dir = fpaths.snapshot_dir;
param_props.snapshot_path = fpaths.snapshot_path;
param_props.figures_path = fpaths.figures_path;

% model path
fpaths.model_path = '/med2/model/';

% temperature and salinity paths
if strcmp(base_grid,'RG')
    fpaths.temp_path = '/med2/gobai_ancillary_data/RG_CLIM/';
    fpaths.sal_path = '/med2/gobai_ancillary_data/RG_CLIM/';
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
    fpaths.param_path = '/fast4/o2/';
    fpaths.data_path = '/data2/GOBAI/o2/';
elseif strcmp(param,'no3')
    fpaths.param_path_temp = '/fast5/no3/';
    fpaths.param_path = '/fast5/no3/';
    fpaths.param_path_o2 = '/fast4/o2/';
    fpaths.data_path = '/data2/GOBAI/no3/';
elseif strcmp(param,'dic')
    fpaths.param_path_temp = '/fast6/dic/';
    fpaths.param_path = '/fast6/dic/';
    fpaths.param_path_o2 = '/fast4/o2/';
    fpaths.param_path_no3 = '/fast5/no3/';
    fpaths.data_path = '/data2/GOBAI/dic/';
else
    error('path_config:UnknownParameter','Unknown parameter: %s',param);
end

param_props.data_path = fpaths.data_path;
param_props.param_path = fpaths.param_path;
param_props.param_path_temp = fpaths.param_path_temp;

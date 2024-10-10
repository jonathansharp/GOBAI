% subsamle_models
%
% DESCRIPTION:
% This function use the combined dataset to subsample models
% for GOBAI validation.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/24/2024

function subsample_models(param,file_date,snap_date,float_file_ext)

%% process parameter name
[param1,param2,param3,~,~,~,~,param_edges] = param_name(param);

%% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load and subsample models
[ESM4,xdim,ydim,zdim,tdim] = load_ESM4('/raid/Model/CMIP6/GFDL-ESM4/',snap_date);

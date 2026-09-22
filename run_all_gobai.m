%% Run all scripts to make all GOBAI data products
t_whole_script=tic; % time entire script

%% Monthly edited configuration parameters
snap_date = 202608;

%% Annually edited configuration parameters
glodap_vrs = 'v3';

%% Other edited configuration parameters
% bounding years
start_year = 1993;
end_year = year(datetime('today'))-1;
% system-specific worker configuration
numWorkers_train = 60;
numWorkers_predict = 60;
numWorkers_cluster = 60;
% float snapshot configuration
snap_download = 1;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
o2_data_modes = {'D' 'A'};
o2_float_file_ext = '_D_A';
no3_data_modes = {'D' 'A'};
no3_float_file_ext = '_D_A';
dic_data_modes = {'D' 'A'};
dic_float_file_ext = '_D_A';
% cluster configuration
num_clusters = 20;
clust_n = 1;
clust_vars = {'temperature_cns' 'salinity_abs' 'pressure'};
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
o2_variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year'};
o2_variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year' 'o2';
o2_variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year' 'o2' 'no3'};
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% data and parameter configuration
data_per_kfold = 0.1; % set data reduction to 10% for k-fold
data_per = 1; % set data reduction to 100% for model training
data_per_osse = 0.2; % set data reduction to 20% for osse
param = 'dic';
param_props = param_config(param);
% base grid
base_grid = 'RFROM';
fpaths = path_config(base_grid,param);
% osse parameters
model_types = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
model_folders = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
realizations = {'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1'};
grid_labels = {'gr' 'gn' 'gn' 'gn' 'gn'};
grid_types = {'regridded' 'native_grid' 'native_grid' 'native_grid' 'native_grid'};
% datasets to include
flt = 1;
gld = 1;
osd = 0;
ctd = 0;

%% plot rfrom animation
plot_rfrom_animation(fpaths,'temp','v2.2_2025','RFROM',...
    start_year,end_year,-1,20,1,'projection','polar','plot_year',2020);
plot_rfrom_animation(fpaths,'sal','v2.2_2025','RFROM',...
    start_year,end_year,-1,20,1,'projection','polar','plot_year',2020);
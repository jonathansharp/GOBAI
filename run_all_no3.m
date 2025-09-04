%% Run all scripts to make GOBAI-NO3
t_whole_script=tic; % time entire script

%% Set configuration parameters
start_year = 1993;
end_year = 2024;
% system-specific worker configuration
numWorkers_train = 20;
numWorkers_predict = 20;
numWorkers_custer = 20;
% float snapshot configuration
snap_download = 1;
snap_date = 202507;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'D'};
float_file_ext = '_D';
% cluster configuration
num_clusters = 15;
clust_vars = {'temperature_cns' 'salinity_abs' 'pressure'};
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year' 'o2'};
% random forest regression configuration
numtrees = 500;
minLeafSize = 10;
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% gradient boosting configuration
numstumps = 500;
numbins = 50;
% data and parameter configuration
data_per_kfold = 0.1; % set data reduction to 10% for k-fold
data_per = 1; % set data reduction to 100% for model training
data_per_osse = 0.15; % set data reduction to 100% for osse
param = 'no3';
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
ctd = 0;

%% load and process data
% acquire data
acquire_snapshot_data(param_props,data_modes,float_file_ext,snap_date,snap_download);
acquire_glodap_data(param_props,glodap_year);
% acquire_wod_ctd_data(param_props,glodap_year,end_year);
% display data
% display_data(param_props,float_file_ext,glodap_year,start_year,snap_date,flt,gld,ctd);
% adjust and combine data
adjust_no3_float_data(float_file_ext,glodap_year,snap_date);
combine_data(param_props,float_file_ext,start_year,glodap_year,snap_date,flt,gld,ctd); % float,glodap,ctd

%% create time-varying clusters and assign data points to them
% form clusters
gmm_clustering(param_props,fpaths,base_grid,start_year,...
    end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
    numWorkers_predict,flt,gld,ctd);
% plot cluster animations
plot_cluster_animation(param_props,fpaths,base_grid,num_clusters,...
    start_year,snap_date,numWorkers_train,flt,gld,ctd);
% plot_probability_animation(base_grid,num_clusters);
% cluster data
assign_data_to_clusters(param_props,base_grid,snap_date,...
    float_file_ext,clust_vars,num_clusters,flt,gld,ctd);
% plot clustered data points
plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,...
    num_clusters,numWorkers_predict,flt,gld,ctd);
% plot_data_over_clusters(param,base_grid,file_date,float_file_ext,...
%    num_clusters,numWorkers_predict);
% develop k-fold evaluation indices
kfold_split_data(param_props,file_date,float_file_ext,...
    num_clusters,num_folds,thresh,flt,gld,ctd);

%% k-fold train models for evaluation statistics
% feed-forward neural networks
train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,flt,gld,ctd,'reduce_data',...
    data_per_kfold,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio,'num_folds',num_folds);

%% train models to create GOBAI product
% feed-forward neural networks
train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,...
    flt,gld,ctd,'reduce_data',data_per,'train_ratio',train_ratio,...
    'val_ratio',val_ratio,'test_ratio',test_ratio);

%% estimate parameter on grid to create GOBAI product
% feed-forward neural networks
predict_gobai('FFNN',param_props,fpaths,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
    end_year,snap_date,flt,gld,ctd,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio);
plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
    file_date,float_file_ext,numWorkers_predict,flt,gld,ctd,'train_ratio',train_ratio,...
    'val_ratio',val_ratio,'test_ratio',test_ratio);

%% run OSSEs
% run_osse(fpaths,model_types,model_folders,realizations,grid_labels,...
%     grid_types,param_props,base_grid,...
%     file_date,snap_date,glodap_year,float_file_ext,start_year,end_year,...
%     num_clusters,variables,clust_vars,train_ratio,val_ratio,test_ratio,...
%     numtrees,minLeafSize,numstumps,numbins,thresh,data_per_osse,...
%     numWorkers_train,numWorkers_predict,flt,gld,ctd);

%% determine uncertainty
% calculate_uncertainty(param_props,base_grid,fpaths,...
%     model_types,num_clusters,numWorkers_predict,file_date,float_file_ext,...
%     glodap_year,train_ratio,val_ratio,test_ratio)
% plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio,'uncer',1);

toc(t_whole_script)

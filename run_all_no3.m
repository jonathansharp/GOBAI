%% Run all scripts to make GOBAI-NO3
t_whole_script=tic; % time entire script

%% Set configuration parameters
start_year = 2004;
end_year = 2022;
% system-specific worker configuration
numWorkers_train = 30;
numWorkers_predict = 30;
% float snapshot configuration
snap_download = 1;
snap_date = 202501;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'D'};
float_file_ext = '_D';
% cluster configuration
num_clusters = 15;
clust_vars = {'temperature_cns' 'salinity_abs' 'pressure'};
glodap_only = false; % EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year'};
% random forest regression configuration
numtrees = 100;
minLeafSize = 10;
% shallow neural network configuration
train_ratio = 0.7;
val_ratio = 0.15;
test_ratio = 0.15;
% gradient boosting configuration
numstumps = 1000;
numbins = 50;
% data and parameter configuration
data_per = 1.0; % set data reduction to 100%
param = 'no3';
param_props = param_config(param);
% base grid
base_grid = 'RFROM';
[model_path,param_path,temp_path,sal_path] = path_config(base_grid,param);

%% load and process data
% acquire data
acquire_snapshot_data(param_props,data_modes,float_file_ext,snap_date,snap_download);
acquire_glodap_data(param_props,glodap_year);
%acquire_wod_ctd_data(param_props,glodap_year);
% display data
display_data(param_props,float_file_ext,glodap_year,snap_date);
% adjust and combine data
adjust_no3_float_data(float_file_ext,glodap_year,snap_date);
combine_data(param_props,float_file_ext,glodap_year,snap_date);

%% create time-varying clusters and assign data points to them
% form clusters
gmm_clustering(param_props,temp_path,sal_path,base_grid,start_year,...
    end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
    numWorkers_predict,param_path);
% plot cluster animations
plot_cluster_animation(param_props,param_path,base_grid,num_clusters,...
    start_year,snap_date,numWorkers_train);
%plot_probability_animation(base_grid,num_clusters);
% cluster data
assign_data_to_clusters(param_props,base_grid,snap_date,...
    float_file_ext,clust_vars,num_clusters);
% plot clustered data points
plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,...
    num_clusters,numWorkers_predict);
%plot_data_over_clusters(param,base_grid,file_date,float_file_ext,...
%    num_clusters,numWorkers_predict);
% develop k-fold evaluation indices
kfold_split_data(param_props,base_grid,file_date,float_file_ext,...
    glodap_only,num_clusters,num_folds,thresh);

%% k-fold train models for evaluation statistics
% feed-forward neural networks
train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio,'num_folds',num_folds);
% random forest regressions
train_gobai('RFR',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'numtrees',numtrees,'minLeafSize',minLeafSize,...
    'num_folds',num_folds);
% gradient-boosting machines
train_gobai('GBM',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'numstumps',numstumps,'numbins',numbins,...
    'num_folds',num_folds);

%% train models to create GOBAI product
% feed-forward neural networks
train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio);
% random forest regressions
train_gobai('RFR',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'numtrees',numtrees,'minLeafSize',minLeafSize);
% gradient-boosting machines
train_gobai('GBM',param_props,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
    data_per,'numstumps',numstumps,'numbins',numbins);

%% estimate parameter on grid to create GOBAI product
% feed-forward neural networks
predict_gobai('FFNN',param_props,param_path,temp_path,sal_path,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
    end_year,snap_date,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio);
plot_gobai_animation(param_props,param_path,base_grid,num_clusters,'FFNN',...
    file_date,float_file_ext,numWorkers_predict,'train_ratio',train_ratio,...
    'val_ratio',val_ratio,'test_ratio',test_ratio);
% random forest regressions
predict_gobai('RFR',param_props,param_path,temp_path,sal_path,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
    end_year,snap_date,'numtrees',numtrees,'minLeafSize',minLeafSize);
plot_gobai_animation(param_props,param_path,base_grid,num_clusters,'RFR',...
    file_date,float_file_ext,numWorkers_predict,'numtrees',numtrees,'minLeafSize',minLeafSize);
% gradient-boosting machines
predict_gobai('GBM',param_props,param_path,temp_path,sal_path,pwd,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
    end_year,snap_date,'numstumps',numstumps,'numbins',numbins);
plot_gobai_animation(param_props,param_path,base_grid,num_clusters,'GBM',...
    file_date,float_file_ext,numWorkers_predict,'numstumps',numstumps,'numbins',numbins);

%% assemble ensemble mean GOBAI
combine_gobai(param_props,temp_path,param_path,base_grid,file_date,float_file_ext,...
    num_clusters,start_year,end_year,snap_date,train_ratio,...
        val_ratio,test_ratio,numtrees,minLeafSize,...
        numstumps,numbins);
plot_gobai_animation(param_props,param_path,base_grid,num_clusters,'AVG',...
    file_date,float_file_ext,numWorkers_predict);

%% run OSSEs
% run_osse(model_path,param_props,file_date,snap_date,float_file_ext,num_clusters,...
%     variables,clust_vars,train_ratio,val_ratio,test_ratio,numtrees,...
%     minLeafSize,numstumps,numbins,thresh,numWorkers_train,numWorkers_predict);

%% determine uncertainty
% calculate_gridding_uncertainty;

toc(t_whole_script)

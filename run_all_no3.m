%% Run all scripts to make GOBAI-NO3
create_config_files; % creates configuration files used by GOBAI code
load_standard_config_files; % loads the 'standard' set of configuration files
load('Config/chinook.mat'); % change 'workers' configuration to chinook
load('Config/base_config_RG.mat'); % change 'base_grid' to Roemmich and Gilson

%% load and process data
% acquire data
acquire_snapshot_data('no3',data_modes,float_file_ext,snap_date,snap_download);
acquire_glodap_data('no3',glodap_year);
% display data
display_data('no3',float_file_ext,file_date,glodap_year);
% adjust and combine data
adjust_no3_float_data(float_file_ext,file_date,glodap_year);
combine_data('no3',float_file_ext,file_date,glodap_year,[0 40],'speed');

%% create time-varying clusters and assign data points to them
% form clusters
gmm_clustering('no3',pwd,base_grid,2004,snap_date,file_date,...
    float_file_ext,clust_vars,num_clusters,numWorkers_predict); %%%%% looked like this worked (12/9/24)
% plot cluster animations
plot_cluster_animation([pwd ''],base_grid,start_year,date_str,num_clusters);
plot_probability_animation(base_grid,num_clusters);
% cluster data
assign_data_to_clusters('no3',base_grid,file_date,snap_date,float_file_ext,clust_vars,num_clusters);
% plot clustered data points
plot_data_by_cluster('no3',base_grid,file_date,float_file_ext,num_clusters);
plot_data_over_clusters('no3',base_grid,file_date,float_file_ext,num_clusters);
% develop k-fold evaluation indices
kfold_split_data('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,num_folds,thresh);

%% k-fold train models for evaluation statistics
% feed-forward neural networks
train_ffnn('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);
% random forest regressions
train_rfr('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);
% gradient-boosting machines
train_gbm('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);

%% train models to create GOBAI product
% feed-forward neural networks
train_ffnn('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);
% random forest regressions
train_rfr('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);
% gradient-boosting machines
train_gbm('no3',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);

%% estimate parameter on grid to create GOBAI product
% feed-forward neural networks
predict_ffnn('no3',pwd,base_grid,file_date,float_file_ext,num_clusters,...
        variables,train_ratio,val_ratio,test_ratio,thresh,...
        numWorkers_predict,2004,snap_date);
gobai_no3_plot_ffnn;
% random forest regressions
gobai_no3_predict_rfr;
gobai_no3_plot_rfr;
% gradient-boosting machines
gobai_no3_predict_gbm;
gobai_no3_plot_gbm;

%% assemble ensemble mean GOBAI
gobai_no3_ensemble;
gobai_no3_plot_ens;

%% run OSSEs
run_osse('no3',file_date,snap_date,float_file_ext,num_clusters,...
    variables,clust_vars,train_ratio,val_ratio,test_ratio,numtrees,...
    minLeafSize,numstumps,numbins,thresh,numWorkers_train,numWorkers_predict);

%% determine uncertainty
% calculate_gridding_uncertainty;

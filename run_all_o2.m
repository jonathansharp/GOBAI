%% Run all scripts to make GOBAI-O2
create_config_files; % creates configuration files used by GOBAI code
load_standard_config_files; % loads the 'standard' set of configuration files
load('Config/chinook.mat'); % change 'workers' configuration to chinook
load('Config/base_config_RG.mat'); % change 'base_grid' to Roemmich and Gilson

%% load and process data
% acquire data
acquire_snapshot_data('o2',data_modes,float_file_ext,snap_date,snap_download);
acquire_glodap_data('o2',glodap_year);
% display data
display_data('o2',float_file_ext,file_date,glodap_year);
% adjust and combine data
adjust_o2_float_data(float_file_ext,file_date,glodap_year);
combine_data('o2',float_file_ext,file_date,glodap_year);

%% create time-varying clusters and assign data points to them
% form clusters
gmm_clustering(pwd,base_grid,2004,snap_date,clust_vars,num_clusters,numWorkers_predict);
% plot cluster animations
plot_cluster_animation([pwd ''],base_grid,start_year,date_str,num_clusters);
plot_probability_animation(base_grid,num_clusters);
% cluster data
assign_data_to_clusters('o2',base_grid,file_date,float_file_ext,clust_vars,num_clusters);
% plot clustered data points
plot_data_by_cluster('o2',base_grid,file_date,float_file_ext,num_clusters);
plot_data_over_clusters('o2',base_grid,file_date,float_file_ext,num_clusters);
% develop k-fold evaluation indices
kfold_split_data('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,num_folds,thresh);

%% k-fold train models for evaluation statistics
% feed-forward neural networks
train_ffnn('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);
% random forest regressions
train_rfr('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);
% gradient-boosting machines
train_gbm('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train,'num_folds',num_folds);

%% train models to create GOBAI product
% feed-forward neural networks
train_ffnn('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);
% random forest regressions
train_rfr('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);
% gradient-boosting machines
train_gbm('o2',base_grid,file_date,float_file_ext,glodap_only,num_clusters,...
    variables,train_ratio,val_ratio,test_ratio,thresh,numWorkers_train);

%% estimate parameter on grid to create GOBAI product
% feed-forward neural networks
predict_ffnn('o2',pwd,base_grid,file_date,float_file_ext,num_clusters,...
        variables,train_ratio,val_ratio,test_ratio,thresh,...
        numWorkers_predict,2004,snap_date);
gobai_o2_plot_ffnn;
% random forest regressions
gobai_o2_predict_rfr;
gobai_o2_plot_rfr;
% gradient-boosting machines
gobai_o2_predict_gbm;
gobai_o2_plot_gbm;

%% assemble ensemble mean GOBAI
gobai_o2_ensemble;
gobai_o2_plot_ens;

%% run OSSEs
run_osse('o2',file_date,snap_date,float_file_ext,num_clusters,...
    variables,clust_vars,train_ratio,val_ratio,test_ratio,numtrees,...
    minLeafSize,numstumps,numbins,thresh,numWorkers_train,numWorkers_predict);

%% determine uncertainty
% calculate_gridding_uncertainty;

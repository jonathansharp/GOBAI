% Run all scripts to make GOBAI-O2
create_config_files;
load_standard_config_files;
load('Config/chinook.mat');

% load data
acquire_snapshot_data('o2',data_modes,float_file_ext,snap_date,snap_download);
acquire_glodap_data('o2',glodap_year);

% display data
display_data('o2',float_file_ext,file_date,glodap_year);

% determine gridding uncertainty
% calculate_gridding_uncertainty;

% adjust and combine data
adjust_o2_float_data(float_file_ext,file_date,glodap_year);
combine_data('o2',float_file_ext,file_date,glodap_year);

% form clusters
gmm_clustering(base_grid,clust_vars,num_clusters,numWorkers_predict);

% plot clusters
plot_cluster_animation(base_grid,num_clusters);
plot_probability_animation(base_grid,num_clusters);

% cluster data
assign_data_to_clusters('o2',base_grid,file_date,float_file_ext,...
    clust_vars,num_clusters);
plot_data_by_cluster('o2',base_grid,file_date,float_file_ext,num_clusters);
plot_data_over_clusters('o2',base_grid,file_date,float_file_ext,num_clusters);
kfold_split_data('o2',base_grid,file_date,float_file_ext,...
    glodap_only,num_clusters,num_folds,thresh);
disp(['data assigned to ' num2str(num_clusters) ' clusters on ' ...
    base_grid ' grid and split into ' num2str(num_folds) ' folds']);

% k-fold train FFNN
dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
    float_file_ext;train_ratio;val_ratio;test_ratio});
kfold_train_ffnn('o2',dir_base,base_grid,file_date,...
    float_file_ext,glodap_only,num_clusters,num_folds,variables,...
    train_ratio,val_ratio,test_ratio,thresh);

gobai_o2_train_ffnn;
gobai_o2_train_rfr;
gobai_o2_train_gbm;
gobai_o2_predict_ffnn;
gobai_o2_plot_ffnn;
gobai_o2_predict_rfr;
gobai_o2_plot_rfr;
gobai_o2_predict_gbm;
gobai_o2_plot_gbm;
gobai_o2_ensemble;

subsample_models('o2',file_date,snap_date,float_file_ext);

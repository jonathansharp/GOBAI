% run_osse
%
% DESCRIPTION:
% This function use the combined dataset to subsample models,
% train test machine learning algorithms, and apply those
% algorithms to model grids to test their performance.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/1/2024

function run_osse(param,file_date,snap_date,float_file_ext,glodap_only,...
    num_clusters,variables,clust_vars,train_ratio,val_ratio,test_ratio,...
    numtrees,minLeafSize,numstumps,numbins,thresh,numWorkers_train,numWorkers_predict)

%% process parameter name
[param1,param2,param3,~,~,~,~,param_edges] = param_name(param);

%% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% define models
model_types = {'GFDL-ESM4'};
model_folders = {'GFDL-ESM4'};

%% loop through each model
for m = 1:length(model_types)

%     %% process models
%     process_cmip_model(model_types{m},['/raid/Model/CMIP6/' model_folders{m} '/'],snap_date,2004,'gr','regridded');
% 
%     %% subsample models
%     subsample_cmip_model(param,param1,all_data,model_types{m},...
%         ['/raid/Model/CMIP6/' model_folders{m} '/'],file_date,...
%         snap_date,float_file_ext,2004,'gr','regridded');
%     
%     %% form clusters
%     gmm_clustering(['/raid/Model/CMIP6/' model_folders{m} '/'],model_types{m},2004,snap_date,clust_vars,num_clusters,numWorkers_predict);
%     % plot_cluster_animation(['/raid/Model/CMIP6/' model_folders{m} '/'],model_types{m},num_clusters,2004,snap_date);
%     
%     %% cluster data
%     assign_data_to_clusters(param,model_types{m},file_date,snap_date,float_file_ext,clust_vars,num_clusters);
%     % plot_data_by_cluster(param,model_types{m},file_date,float_file_ext,num_clusters,numWorkers_train);
%     % plot_data_over_clusters(param,model_types{m},file_date,float_file_ext,num_clusters,numWorkers_train);
%     
%     %% train algorithms
%     data_per = 0.02; % fraction for reducing training data volume
%     % train feed forward neural networks
%     train_ffnn(param,model_types{m},file_date,float_file_ext,...
%         num_clusters,variables,train_ratio,val_ratio,test_ratio,thresh,...
%         numWorkers_train,'reduce_data',data_per);
%     % train random forest regressions
%     train_rfr(param,model_types{m},file_date,float_file_ext,...
%         num_clusters,variables,numtrees,minLeafSize,thresh,...
%         numWorkers_train,'reduce_data',data_per);
%     % train gradient-boosting machines
%     train_gbm(param,model_types{m},file_date,float_file_ext,...
%         num_clusters,variables,numstumps,numbins,thresh,...
%         numWorkers_train,'reduce_data',data_per);
% 
%     %% predict on grid
%     % predict on grid using feed forward neural networks
%     predict_gobai('FFNN',param,['/raid/Model/CMIP6/' model_folders{m} '/'],...
%         model_types{m},file_date,float_file_ext,num_clusters,variables,thresh,...
%         numWorkers_predict,2004,snap_date,'train_ratio',train_ratio,...
%         'val_ratio',val_ratio,'test_ratio',test_ratio);
%     % predict on grid using random forest regressions
%     predict_gobai('RFR',param,['/raid/Model/CMIP6/' model_folders{m} '/'],...
%         model_types{m},file_date,float_file_ext,num_clusters,...
%         variables,thresh,numWorkers_predict,2004,snap_date,...
%         'numtrees',numtrees,'minLeafSize',minLeafSize);
%     % predict on grid using gradient-boosting machines
%     predict_gobai('GBM',param,['/raid/Model/CMIP6/' model_folders{m} '/'],...
%         model_types{m},file_date,float_file_ext,num_clusters,...
%         variables,thresh,numWorkers_predict,2004,snap_date,...
%         'numstumps',numstumps,'numbins',numbins);

    %% average grids across all algorithms
    combine_gobai(param,model_types{m},file_date,float_file_ext,...
        num_clusters,numWorkers_predict,2004,snap_date,train_ratio,...
    val_ratio,test_ratio,numtrees,minLeafSize,numstumps,numbins);

end

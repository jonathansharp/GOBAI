% run_osse
%
% DESCRIPTION:
% This function use the combined dataset to subsample models,
% train test machine learning algorithms, and apply those
% algorithms to model grids to test their performance.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 4/15/2025

function run_osse(model_path,param_props,base_grid,param_path,file_date,snap_date,glodap_year,float_file_ext,...
    start_year,end_year,num_clusters,variables,clust_vars,train_ratio,val_ratio,test_ratio,...
    numtrees,minLeafSize,numstumps,numbins,thresh,numWorkers_train,numWorkers_predict)

%% load combined data
load([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% fraction for reducing training data volume
data_per = 1; % 100%

%% define models
model_types = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
model_folders = {'GFDL-ESM4' 'CanESM5' 'IPSL-CM6A-LR' 'ACCESS-ESM1-5' 'MPI-ESM1-2-LR'};
realizations = {'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1' 'r1i1p1f1'};
grid_labels = {'gr' 'gn' 'gn' 'gn' 'gn'};
grid_types = {'regridded' 'native_grid' 'native_grid' 'native_grid' 'native_grid'};

%% loop through each model
for m = 1:length(model_types)

    %% process models
    process_cmip_model(model_types{m},[model_path model_folders{m} '/'],...
        snap_date,2004,realizations{m},grid_labels{m},grid_types{m});

    %% subsample models
    subsample_cmip_model(param_props,all_data,model_types{m},...
        [model_path model_folders{m} '/'],file_date,...
        snap_date,float_file_ext,2004,realizations{m});

    %% form clusters
    gmm_clustering(param_props,[model_path model_folders{m} '/'],...
        [model_path model_folders{m} '/'],model_types{m},...
        start_year,end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
        numWorkers_predict,[model_path model_folders{m} '/'],'rlz',realizations{m});

    %% cluster data
    assign_data_to_clusters(param_props,model_types{m},snap_date,float_file_ext,clust_vars,num_clusters);
    plot_data_by_cluster(param_props,model_types{m},file_date,float_file_ext,num_clusters,numWorkers_train);
    % plot_data_over_clusters(param,model_types{m},file_date,float_file_ext,num_clusters,numWorkers_train);

    %% train algorithms
    % train feed forward neural networks
    train_gobai('FFNN',param_props,model_types{m},file_date,float_file_ext,...
        num_clusters,variables,thresh,numWorkers_train,snap_date,'reduce_data',...
        data_per,'train_ratio',train_ratio,'val_ratio',val_ratio,...
        'test_ratio',test_ratio);

    %% predict on grid
    % predict on grid using feed forward neural networks
    predict_gobai('FFNN',param_props,[model_path model_folders{m} '/'],...
        [model_path model_folders{m} '/'],[model_path model_folders{m} '/'],...
        model_types{m},file_date,float_file_ext,num_clusters,...
        variables,thresh,numWorkers_predict,clust_vars,start_year,end_year,...
        snap_date,'train_ratio',train_ratio,'val_ratio',val_ratio,...
        'test_ratio',test_ratio,'rlz',realizations{m});

    %% average grids across all algorithms
    % combine_gobai(param_props,[model_path model_folders{m} '/'],...
    %     [model_path model_folders{m} '/'],...
    %     model_types{m},file_date,float_file_ext,...
    %     num_clusters_1,num_clusters_2,num_clusters_3,start_year,end_year,snap_date,train_ratio,...
    %     val_ratio,test_ratio,numtrees,minLeafSize,...
    %     numstumps,numbins,'rlz',realizations{m});

    %% compare reconstructed model grid to original
    compare_osse(param_props,[model_path model_folders{m} '/'],...
        model_types{m},file_date,float_file_ext,...
        num_clusters,2004,snap_date,train_ratio,...
        val_ratio,test_ratio,realizations{m});

end

% compare all osse experiments
compare_all_osses(param_props,model_path,model_types,realizations,...
    num_clusters,file_date,float_file_ext,train_ratio,val_ratio,test_ratio);

% calculate uncertainty for gobai
calculate_uncertainty(param_props,base_grid,param_path,model_path,model_types,realizations,...
    num_clusters,numWorkers_predict,file_date,float_file_ext,glodap_year,train_ratio,val_ratio,test_ratio)


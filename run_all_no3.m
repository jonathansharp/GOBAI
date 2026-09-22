%% Run all scripts to make GOBAI-NO3
t_whole_script=tic; % time entire script

%% Set configuration parameters
start_year = 1993;
end_year = 2025;
% end_year = year(datetime('today'));
% system-specific worker configuration
numWorkers_train = 24;
numWorkers_predict = 24;
numWorkers_cluster = 24;
% float snapshot configuration
snap_download = 1;
snap_date = 202606;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_vrs = 'v3';
data_modes = {'D' 'A'};
float_file_ext = '_D_A';
% cluster configuration
num_clusters = 20;
clust_n = 1;
clust_vars = {'temperature_cns' 'salinity_abs' 'pressure'};
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year' 'o2'};
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% data and parameter configuration
data_per_kfold = 0.1; % set data reduction to 10% for k-fold
data_per = 1; % set data reduction to 100% for model training
data_per_osse = 0.2; % set data reduction to 20% for osse
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
osd = 1;
ctd = 0;

%% load and process data
% % acquire data
% acquire_snapshot_data(param_props,data_modes,float_file_ext,snap_date,snap_download);
% acquire_glodap_data(param_props,glodap_vrs,start_year);
% if osd == 1 | ctd == 1; acquire_wod_data(param_props,...
%     glodap_vrs,start_year,osd,ctd,numWorkers_predict); end
% % adjust data
% adjust_no3_data(float_file_ext,glodap_vrs,snap_date,flt,gld,osd,ctd);
% % display data
% display_data(param_props,float_file_ext,glodap_vrs,...
%     start_year,snap_date,flt,gld,osd,ctd);
% % combine data
% combine_data(param_props,float_file_ext,...
%     glodap_vrs,snap_date,flt,gld,osd,ctd); % float,glodap,ctd
% 
% %% plot histogram of data
% plot_data_hist(param_props,file_date,float_file_ext,...
%     flt,gld,osd,ctd,start_year);
% 
% %% determine ideal number of clusters
% num_clusters = [10;12;14;16;18;20;22;24;26;28;30];
% RMSE = nan(size(num_clusters));
% MAE = nan(size(num_clusters));
% R2 = nan(size(num_clusters));
% kfold_table = table(num_clusters,RMSE,MAE,R2);
% 
% for clust_n = 1:length(num_clusters)
% 
% %% create time-varying clusters and assign data points to them
% % form clusters
% gmm_clustering(param_props,fpaths,base_grid,start_year,...
%     end_year,snap_date,float_file_ext,clust_vars,num_clusters(clust_n),...
%     numWorkers_predict,flt,gld,osd,ctd);
% % % plot cluster animations
% % plot_cluster_animation(param_props,fpaths,base_grid,num_clusters(clust_n),...
% %     start_year,snap_date,numWorkers_train,flt,gld,ctd);
% % plot_probability_animation(base_grid,num_clusters);
% % cluster data
% assign_data_to_clusters(param_props,base_grid,snap_date,...
%     float_file_ext,clust_vars,num_clusters(clust_n),flt,gld,osd,ctd);
% % plot clustered data points
% plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,...
%     num_clusters(clust_n),numWorkers_predict,flt,gld,osd,ctd);
% % plot_data_over_clusters(param,base_grid,file_date,float_file_ext,...
% %    num_clusters,numWorkers_predict);
% % develop k-fold evaluation indices
% kfold_split_data(param_props,file_date,float_file_ext,...
%     num_clusters(clust_n),num_folds,thresh,flt,gld,osd,ctd);
% 
% %% k-fold train models for evaluation statistics
% % feed-forward neural networks
% train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
%     num_clusters(clust_n),variables,thresh,numWorkers_train,snap_date,...
%     flt,gld,osd,ctd,'reduce_data',data_per_kfold,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio,'num_folds',num_folds);
% 
% %% assemble table of K-fold error statistics
% kfold_dir = [param_props.dir_name '/KFold/FFNN/c' ...
%     num2str(num_clusters(clust_n)) '_' file_date float_file_ext];
% kfold_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' ...
%     num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
% load([kfold_dir '/' kfold_name],'alg_rmse','alg_med_abs_err','alg_r2');
% kfold_table(clust_n,'RMSE') = {alg_rmse};
% kfold_table(clust_n,'MAE') = {alg_med_abs_err};
% kfold_table(clust_n,'R2') = {alg_r2};
% clear alg_rmse alg_med_abs_err alg_r2}
% 
% end
% 
% % save k-fold table
% save([param_props.dir_name '/KFold/k_fold_stats_' ...
%     char(datetime('today','Format','dd-MMM-uuuu'))],'kfold_table');
% clear kfold_table

%% train models to create GOBAI product
% % feed-forward neural networks
% train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_train,snap_date,...
%     flt,gld,osd,ctd,'reduce_data',data_per,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio);

%% estimate parameter on grid to create GOBAI product
% % feed-forward neural networks
% predict_gobai('FFNN',param_props,fpaths,base_grid,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
%     end_year,snap_date,flt,gld,osd,ctd,'train_ratio',train_ratio,'val_ratio',val_ratio,...
%     'test_ratio',test_ratio);
% % plot gobai animations
% plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,flt,gld,osd,ctd,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio);
% plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,flt,gld,ctd,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio,'anom',1);
% 
% % plot gobai timeseries
% gobai_global_timeseries(base_grid,param_props,fpaths,num_clusters,file_date,...
%     float_file_ext,train_ratio,val_ratio,test_ratio,flt,gld,ctd,start_year,end_year)

%% determine ideal number of clusters
num_clusters = [10;12;14;16;18;20;22;24;26;28;30];
RMSE = nan(size(num_clusters));
MAE = nan(size(num_clusters));
R2 = nan(size(num_clusters));
kfold_table = table(num_clusters,RMSE,MAE,R2);
for clust_n = 1:length(num_clusters)

%% run OSSEs
run_osse(fpaths,model_types,model_folders,realizations,grid_labels,...
    grid_types,param_props,base_grid,file_date,snap_date,glodap_vrs,...
    float_file_ext,start_year,end_year,num_clusters(clust_n),variables,...
    clust_vars,train_ratio,val_ratio,test_ratio,thresh,data_per_osse,...
    numWorkers_train,numWorkers_predict,flt,gld,osd,ctd);

% compare_two_osse_ensembles(param_props,fpaths,model_types,realizations,...
%     num_clusters,file_date,float_file_ext,train_ratio,val_ratio,test_ratio,...
%     start_year,snap_date,'g','fg');

end

%% assemble table of K-fold error statistics
% kfold_dir = [param_props.dir_name '/KFold/FFNN/c' ...
%     num2str(num_clusters(clust_n)) '_' file_date float_file_ext];
% kfold_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' ...
%     num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
% load([kfold_dir '/' kfold_name],'alg_rmse','alg_med_abs_err','alg_r2');
% kfold_table(clust_n,'RMSE') = {alg_rmse};
% kfold_table(clust_n,'MAE') = {alg_med_abs_err};
% kfold_table(clust_n,'R2') = {alg_r2};
% clear alg_rmse alg_med_abs_err alg_r2}

% compare_two_osse_ensembles(param_props,fpaths,model_types,realizations,...
%     num_clusters,file_date,float_file_ext,train_ratio,val_ratio,test_ratio,...
%     start_year,snap_date,'g','fg');

%% determine uncertainty
% calculate_uncertainty(param_props,base_grid,fpaths,...
%     model_types,num_clusters,numWorkers_predict,file_date,float_file_ext,...
%     glodap_vrs,train_ratio,val_ratio,test_ratio)
% plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio,'uncer',1);

toc(t_whole_script)



%% Set configuration parameters
start_year = 1993;
end_year = 2025;
% end_year = year(datetime('today'));
% system-specific worker configuration
numWorkers_train = 24;
numWorkers_predict = 24;
% float snapshot configuration
snap_download = 1;
snap_date = 202602;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'D' 'A'};
float_file_ext = '_D_A';
% cluster configuration
num_clusters = 15;
clust_n = 1;
clust_vars = {'temperature_cns' 'salinity_abs' 'sigma'};
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year'};
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% data and parameter configuration
data_per_kfold = 0.2; % set data reduction to 10% for k-fold
data_per = 1; % set data reduction to 100% for model training
data_per_osse = 0.01; % set data reduction to 10% for osse
param = 'o2';
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

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% load combined data
load([param_props.dir_name '/Data/processed_all_' param_props.file_name ...
    '_data_' float_ext glodap_ext ctd_ext '_' file_date float_file_ext ...
    '.mat'],'all_data','file_date');

%% just for IPSL-model
m = 3;

%% process models
process_cmip_model(model_types{m},[fpaths.model_path model_folders{m} '/'],...
    snap_date,start_year,realizations{m},grid_labels{m},grid_types{m});

%% extract 2025 coverage
all_data_2025_cov = struct();
idx_2025 = all_data.year == 2025;
vars = fieldnames(all_data);
for v = 1:length(vars)
    if ndims(all_data.(vars{v})) <= 2
        all_data_2025_cov.(vars{v}) = all_data.(vars{v})(idx_2025);
    end
end

%% replicate 2025 coverage
all_data_2025_cov.year = repelem((1993:2025),length(all_data_2025_cov.year))';
vars = fieldnames(all_data_2025_cov);
for v = 1:length(vars)
    if ndims(all_data_2025_cov.(vars{v})) <= 2 & ~strcmp((vars{v}),'year')
        all_data_2025_cov.(vars{v}) = repmat(all_data_2025_cov.(vars{v}),...
            max(all_data_2025_cov.year)-min(all_data_2025_cov.year)+1,1);
    end
end
all_data_2025_cov.time = datenum(all_data_2025_cov.year,1,all_data_2025_cov.day);

%% subsample models
% normal coverage
subsample_cmip_model(param_props,all_data,model_types{m},...
    [fpaths.model_path model_folders{m} '/'],file_date,...
    snap_date,float_file_ext,start_year,realizations{m},...
    float_ext,glodap_ext,ctd_ext);
% 2025 coverage
subsample_cmip_model(param_props,all_data_2025_cov,model_types{m},...
    [fpaths.model_path model_folders{m} '/'],file_date,...
    snap_date,float_file_ext,start_year,realizations{m},...
    float_ext,glodap_ext,ctd_ext,'coverage',2025);

%% form clusters
% normal coverage
gmm_clustering(param_props,fpaths,model_types{m},...
    start_year,end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
    numWorkers_predict,flt,gld,ctd,'rlz',realizations{m});
% 2025 coverage
gmm_clustering(param_props,fpaths,model_types{m},...
    start_year,end_year,snap_date,float_file_ext,clust_vars,num_clusters,...
    numWorkers_predict,flt,gld,ctd,'rlz',realizations{m},'coverage',2025);

%% cluster data
% normal coverage
assign_data_to_clusters(param_props,model_types{m},snap_date,...
    float_file_ext,clust_vars,num_clusters,flt,gld,ctd);
plot_data_by_cluster(param_props,model_types{m},file_date,...
    float_file_ext,num_clusters,numWorkers_train,flt,gld,ctd);
% 2025 coverage
assign_data_to_clusters(param_props,model_types{m},snap_date,...
    float_file_ext,clust_vars,num_clusters,flt,gld,ctd,...
    'coverage',2025);
plot_data_by_cluster(param_props,model_types{m},file_date,...
    float_file_ext,num_clusters,numWorkers_train,flt,gld,ctd,...
    'coverage',2025);

%% train algorithms
% normal coverage
train_gobai('FFNN',param_props,model_types{m},file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,...
    flt,gld,ctd,'reduce_data',data_per_osse,'train_ratio',...
    train_ratio,'val_ratio',val_ratio,'test_ratio',test_ratio);
% 2025 coverage
train_gobai('FFNN',param_props,model_types{m},file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_train,snap_date,...
    flt,gld,ctd,'reduce_data',data_per_osse,'train_ratio',...
    train_ratio,'val_ratio',val_ratio,'test_ratio',test_ratio,...
    'coverage',2025);

%% predict on grid
% normal coverage
predict_gobai('FFNN',param_props,fpaths,...
    model_types{m},file_date,float_file_ext,num_clusters,...
    variables,thresh,numWorkers_predict,clust_vars,start_year,end_year,...
    snap_date,flt,gld,ctd,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio,'rlz',realizations{m});
% 2025 coverage
predict_gobai('FFNN',param_props,fpaths,...
    model_types{m},file_date,float_file_ext,num_clusters,...
    variables,thresh,numWorkers_predict,clust_vars,start_year,end_year,...
    snap_date,flt,gld,ctd,'train_ratio',train_ratio,'val_ratio',val_ratio,...
    'test_ratio',test_ratio,'rlz',realizations{m},'coverage',2025);

%% compare reconstructed model grid to original
% normal coverage
compare_osse(param_props,fpaths,...
    model_types{m},file_date,float_file_ext,...
    num_clusters,start_year,snap_date,train_ratio,...
    val_ratio,test_ratio,realizations{m},...
    float_ext,glodap_ext,ctd_ext,numWorkers_predict,...
    'coverage',2025);
% 2025 coverage
compare_osse(param_props,fpaths,...
    model_types{m},file_date,float_file_ext,...
    num_clusters,start_year,snap_date,train_ratio,...
    val_ratio,test_ratio,realizations{m},...
    float_ext,glodap_ext,ctd_ext,numWorkers_predict,...
    'coverage',2025);
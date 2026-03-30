%% Run all scripts to make GOBAI-DIC
t_whole_script=tic; % time entire script

%% Set configuration parameters
start_year = 1993;
end_year = 2025;
% system-specific worker configuration
numWorkers_train = 24;
numWorkers_predict = 24;
numWorkers_cluster = 24;
% float snapshot configuration
snap_download = 1;
snap_date = 202603;
file_date = datestr(datenum(floor(snap_date/1e2),...
    mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'D'};
float_file_ext = '_D';
% cluster configuration
num_clusters = 15;
clust_n = 1;
clust_vars = {'temperature_cns' 'salinity_abs' 'sigma'};
thresh = 0.05;
num_folds = 5;
% algorithm training configuration
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year' 'o2' 'no3'};
% shallow neural network configuration
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
% data and parameter configuration
data_per_kfold = 0.2; % set data reduction to 20% for k-fold
data_per = 1; % set data reduction to 100%
data_per_osse = 1; % set data reduction to 20% for osse
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
ctd = 0;

%% plot rfrom animation
% plot_rfrom_temp_animation(fpaths,'v2.2','RFROM',start_year,end_year)
% plot_rfrom_sal_animation(fpaths,'v2.2','RFROM',start_year,end_year)

%% load and process data
% % acquire data
% acquire_snapshot_data(param_props,data_modes,float_file_ext,snap_date,0);
% acquire_glodap_data(param_props,glodap_year,start_year);
% % display data
% display_data(param_props,float_file_ext,glodap_year,start_year,snap_date,flt,gld,ctd);
% % adjust and combine data
% if flt == 1; adjust_dic_float_data(float_file_ext,glodap_year,snap_date); end
% combine_data(param_props,float_file_ext,start_year,glodap_year,snap_date,flt,gld,ctd);

%% plot histogram of data
% plot_data_hist(param_props,file_date,float_file_ext,...
%     flt,gld,ctd,start_year,end_year);

%% determine ideal number of clusters
% num_clusters = [20 22 25 27 30];
% for clust_n = 1:length(num_clusters)

%% create time-varying clusters and assign data points to them
% form clusters
% gmm_clustering(param_props,fpaths,base_grid,start_year,...
%     end_year,snap_date,float_file_ext,clust_vars,num_clusters(clust_n),...
%     numWorkers_predict,flt,gld,ctd);
% % % plot cluster animations
% % plot_cluster_animation(param_props,fpaths,base_grid,num_clusters(clust_n),...
% %     start_year,snap_date,numWorkers_train,flt,gld,ctd);
% % plot_probability_animation(base_grid,num_clusters);
% % cluster data
% assign_data_to_clusters(param_props,base_grid,snap_date,...
%     float_file_ext,clust_vars,num_clusters(clust_n),flt,gld,ctd);
% % plot clustered data points
% plot_data_by_cluster(param_props,base_grid,file_date,float_file_ext,...
%     num_clusters(clust_n),numWorkers_predict,flt,gld,ctd);
% % plot_data_over_clusters(param,base_grid,file_date,float_file_ext,...
% %    num_clusters,numWorkers_predict);
% % develop k-fold evaluation indices
% kfold_split_data(param_props,file_date,float_file_ext,...
%     num_clusters(clust_n),num_folds,thresh,flt,gld,ctd);

%% k-fold train models for evaluation statistics
% feed-forward neural networks
train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
    num_clusters(clust_n),variables,thresh,numWorkers_train,snap_date,...
    flt,gld,ctd,'reduce_data',data_per_kfold,'train_ratio',train_ratio,...
    'val_ratio',val_ratio,'test_ratio',test_ratio,'num_folds',num_folds);

% end

%% train models to create GOBAI product
% % feed-forward neural networks
% train_gobai('FFNN',param_props,base_grid,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_train,snap_date,...
%     flt,gld,ctd,'reduce_data',data_per,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio);

%% estimate parameter on grid to create GOBAI product
% % feed-forward neural networks
% predict_gobai('FFNN',param_props,fpaths,base_grid,file_date,float_file_ext,...
%     num_clusters,variables,thresh,numWorkers_predict,clust_vars,start_year,...
%     end_year,snap_date,flt,gld,ctd,'train_ratio',train_ratio,'val_ratio',val_ratio,...
%     'test_ratio',test_ratio);
% % plot gobai animations
% plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
%     file_date,float_file_ext,numWorkers_predict,flt,gld,ctd,'train_ratio',train_ratio,...
%     'val_ratio',val_ratio,'test_ratio',test_ratio);
% % plot_gobai_animation(param_props,fpaths,base_grid,num_clusters,'FFNN',...
% %     file_date,float_file_ext,numWorkers_predict,flt,gld,ctd,'train_ratio',train_ratio,...
% %     'val_ratio',val_ratio,'test_ratio',test_ratio,'anom',1);

% % plot gobai timeseries
% gobai_global_timeseries(base_grid,param_props,fpaths,num_clusters,file_date,...
%     float_file_ext,train_ratio,val_ratio,test_ratio,flt,gld,ctd,start_year,end_year)

%% run OSSEs
run_osse(fpaths,model_types,model_folders,realizations,grid_labels,...
    grid_types,param_props,base_grid,...
    file_date,snap_date,glodap_year,float_file_ext,start_year,end_year,...
    num_clusters(clust_n),variables,clust_vars,train_ratio,val_ratio,test_ratio,...
    numtrees,minLeafSize,numstumps,numbins,thresh,data_per_osse,...
    numWorkers_train,numWorkers_predict,flt,gld,ctd);

%% determine uncertainty
% calculate_gridding_uncertainty;

%% evaluate timeseries
% lon_fig = 200;
% lat_fig = 40;
% pres_fig = 10;
% 
% filename = [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
%             '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
%             '_val' num2str(100*test_ratio) '_test' num2str(100*val_ratio) ...
%             '/gobai-' param_props.file_name '.nc'];
% lon = ncread(filename,'lon'); [~,lon_idx] = min(abs(lon-lon_fig));
% lat = ncread(filename,'lat'); [~,lat_idx] = min(abs(lat-lat_fig));
% pres = ncread(filename,'pres'); [~,pres_idx] = min(abs(pres-pres_fig));
% 
% % global mean dic
% for t = 1:992
%     dic_temp = ncread(filename,'dic',[1 1 1 t],[Inf Inf Inf 1]);
%     dic(t) = mean(dic_temp(:),'omitnan');
% end
% 
% % dic = squeeze(ncread(filename,'dic',[lon_idx lat_idx pres_idx 1],[1 1 1 Inf]));
% time = ncread(filename,'time');
% 
% figure; hold on; set(gcf,'position',[100 100 1000 200]);
% dic_m = movmean(dic,52);
% mdl=fitlm(time,dic_m);
% dic_fit = (mdl.Coefficients{2,1}.*time + mdl.Coefficients{1,1})';
% ylabel('DIC Anomaly (\mumol kg^{-1})');
% clr=cmocean('amplitude',1);
% plot(time(53:end-52),dic_m(53:end-52)-dic_fit(53:end-52),'LineWidth',3,'Color',clr);
% % plot(time,dic-mean(dic),'LineWidth',1,'LineStyle','-','Color','#0072BD');
% plot([time(53) time(end-52)],[0 0],'LineWidth',1,'LineStyle',':','Color','k');
% datetick('x');

%% end timing
toc(t_whole_script)

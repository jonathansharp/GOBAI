% kfold_avg_all
%
% DESCRIPTION:
% This function takes an average of all k-fold evaluation results.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 4/24/2025

function kfold_avg_all(param_props,base_grid,float_file_ext,num_clusters,...
    snap_date,train_ratio,val_ratio,test_ratio,numtrees,minLeafSize,...
    numstumps,numbins,varargin)

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% create directory names
% FFNN
kfold_ffnn_dir = ...
    ['KFold/FFNN/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_ffnn_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*test_ratio)];
% GBM
kfold_gbm_dir = ...
    ['KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_gbm_name = ['GBM_output_tr' num2str(numstumps) '_bin' num2str(numbins)];
% RFR
kfold_rfr_dir = ...
    ['KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_rfr_name = ['RFR_output_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
% AVG
kfold_ens_dir = ...
    ['KFold/AVG/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_ens_name = ['AVG_output' ...
    'FFNN_train' num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' num2str(100*val_ratio) ...
    'GBM_tr' num2str(numstumps) 'RFR_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
% Figures
fig_dir = [param_props.dir_name '/Figures/KFold/AVG/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name = 'k_fold_comparison.png';
fig_name_2 = 'k_fold_spatial_comparison.png';

%% load k-fold data
% FFNN
ffnn_output = load([param_props.dir_name '/' kfold_ffnn_dir '/' kfold_ffnn_name '.mat']);
% GBM
gbm_output = load([param_props.dir_name '/' kfold_gbm_dir '/' kfold_gbm_name '.mat']);
% RFR
rfr_output = load([param_props.dir_name '/' kfold_rfr_dir '/' kfold_rfr_name '.mat']);

%% average across models
ens_output.(['k_fold_test_' param_props.file_name]) = ...
    mean([rfr_output.alg_output.(['k_fold_test_' param_props.file_name]) ...
          ffnn_output.alg_output.(['k_fold_test_' param_props.file_name]) ...
          gbm_output.alg_output.(['k_fold_test_' param_props.file_name])],2);
% clean up
clear rfr_output ffnn_output gbm_output
% compare k-fold output to data
ens_output.k_fold_delta = ...
    ens_output.(['k_fold_test_' param_props.file_name]) - all_data.(param_props.file_name);
% calculate error stats
ens_mean_err = mean(ens_output.k_fold_delta,'omitnan');
ens_med_err = median(ens_output.k_fold_delta,'omitnan');
ens_rmse = sqrt(mean(ens_output.k_fold_delta.^2,'omitnan'));
ens_med_abs_err = median(abs(ens_output.k_fold_delta),'omitnan');
% print error stats
fprintf(['Mean Error = ' num2str(ens_mean_err) ' ' param_props.units '\n']);
fprintf(['Median Error = ' num2str(ens_med_err) ' ' param_props.units '\n']);
fprintf(['RMSE = ' num2str(ens_rmse) ' ' param_props.units '\n']);
fprintf(['Median Abs. Error = ' num2str(ens_med_abs_err) ' ' param_props.units '\n']);
% save predicted data
if ~isfolder([pwd '/' kfold_ens_dir]); mkdir(kfold_ens_dir); end
save([kfold_ens_dir '/' kfold_ens_name],...
    'ens_output','ens_rmse','ens_med_err','ens_mean_err');

%% plot histogram of errors
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.(param_props.file_name) ...
    ens_output.(['k_fold_test_' param_props.file_name])],...
    'Edges',{param_props.edges param_props.edges});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([param_props.edges(1) param_props.edges(end)],[param_props.edges(1) param_props.edges(end)],'k--');
set(h,'EdgeColor','none');
xlim([param_props.edges(1) param_props.edges(end)]);
ylim([param_props.edges(1) param_props.edges(end)]);
xlabel(['Measured ' param_props.label ' ' param_props.units]);
ylabel(['AVG ' param_props.label ' ' param_props.units]);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
clim([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(param_props.edges(1)+(2/5)*(param_props.edges(end)-param_props.edges(1)),...
    param_props.edges(1)+(1/10)*(param_props.edges(end)-param_props.edges(1)),...
    ['RMSE = ' num2str(round(ens_rmse,1)) ' ' param_props.units],'fontsize',12);
if ~isfolder(fig_dir); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name]);
% clean up
clear counts bin_centers h p myColorMap
close

%% plot gridded errors
% determine bin number of each test data point on 1 degree grid
lon_edges = -180:180; lon = -179.5:179.5;
lat_edges = -90:90; lat = -89.5:89.5;
[~,~,Xnum] = histcounts(all_data.longitude,lon_edges);
[~,~,Ynum] = histcounts(all_data.latitude,lat_edges);
% accumulate 3D grid of test data point errors
subs = [Xnum, Ynum];
idx_subs = any(subs==0,2);
sz = [length(lon),length(lat)];
ens_output.k_fold_delta_spatial = accumarray(subs(~idx_subs,:),...
    abs(ens_output.k_fold_delta(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure('visible','off'); hold on
worldmap([-90 90],[20 380]);
setm(gca,'mapprojection','robinson');
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(lat,[lon lon(end)+1],[ens_output.k_fold_delta_spatial ...
    ens_output.k_fold_delta_spatial(:,end)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
clim([0 (2/50)*(param_props.edges(end)-param_props.edges(1))]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Absolute \Delta' param_props.label];
c.FontSize = 22;
c.TickLength = 0;
mlabel off; plabel off;
if ~isfolder(fig_dir); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_2]);
% clean up
clear land cmap c
close

%% clean up
clear

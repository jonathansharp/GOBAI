% k_fold_ensemble
%
% DESCRIPTION:
% This function uses validation versions of machine learning models to
% evaluate a reserved subset of data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/2/2024

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% create directory names
% FFNN
kfold_ffnn_dir = ...
    ['KFold/FFNN/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_ffnn_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
% GBM
kfold_gbm_dir = ...
    ['KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_gbm_name = ['GBM_output_tr' num2str(numstumps)];
% RFR
kfold_rfr_dir = ...
    ['KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_rfr_name = ['RFR_output_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
% FFNN
kfold_ens_dir = ...
    ['KFold/ENS/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_ens_name = ['ENS_output' ...
    'FFNN_train' num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' num2str(100*val_ratio) ...
    'GBM_tr' num2str(numstumps) 'RFR_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
% Figures
fig_dir = ['Figures/KFold/ENS/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name = ['k_fold_comparison.png'];

%% load k-fold data
% FFNN
load([kfold_ffnn_dir '/' kfold_ffnn_name],'ffnn_output','ffnn_rmse',...
    'ffnn_med_err','ffnn_mean_err');
% GBM
load([kfold_gbm_dir '/' kfold_gbm_name],'gbm_output','gbm_rmse',...
    'gbm_med_err','gbm_mean_err');
% RFR
load([kfold_rfr_dir '/' kfold_rfr_name],'rfr_output','rfr_rmse',...
    'rfr_med_err','rfr_mean_err');

%% average across models
ens_output.k_fold_test_oxygen = ...
    mean([rfr_output.k_fold_test_oxygen ffnn_output.k_fold_test_oxygen ...
    gbm_output.k_fold_test_oxygen],2);
% clean up
clear rfr_output ffnn_output gbm_output
% compare k-fold output to data
ens_output.k_fold_delta = ens_output.k_fold_test_oxygen - all_data.oxygen;
% calculate error stats
ens_mean_err = mean(ens_output.k_fold_delta);
ens_med_err = median(ens_output.k_fold_delta);
ens_rmse = sqrt(mean(ens_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_ens_dir]); mkdir(kfold_ens_dir); end
save([kfold_ens_dir '/' kfold_ens_name],...
    'ens_output','ens_rmse','ens_med_err','ens_mean_err');

%% plot histogram of errors
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.oxygen ens_output.k_fold_test_oxygen],...
    'Edges',{0:5:500 0:5:500});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([0 500],[0 500],'k--');
set(h,'EdgeColor','none');
xlim([0 500]); ylim([0 500]);
xlabel('Measured Oxygen (\mumol kg^{-1})');
ylabel('ENS Oxygen (\mumol kg^{-1})');
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(300,50,['RMSE = ' num2str(round(ens_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
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
figure; hold on
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
caxis([0 20]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Absolute \Delta[O_{2}]'];
c.FontSize = 22;
c.TickLength = 0;
mlabel off; plabel off;
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_2]);
% clean up
clear land cmap c
close

%% clean up
clear

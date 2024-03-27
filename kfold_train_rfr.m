% k_fold_fit_test_models
%
% DESCRIPTION:
% This function uses a subset of the combined dataset to train validation
% versions of machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/15/2024

function kfold_train_rfr(param,dir_base,base_grid,file_date,...
    float_file_ext,glodap_only,num_clusters,num_folds,variables,...
    numtrees,minLeafSize,thresh)

%% process parameter name
[param1,param2,param3,~,~,~,~,param_edges] = param_name(param);

%% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load([param1 '/Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% load data cluster indices
load([param1 '/Data/k_fold_data_indices_'  base_grid '_' num2str(num_clusters) ...
    '_' num2str(num_folds) '_' file_date float_file_ext '.mat'],...
    'num_folds','train_idx','test_idx');

%% remove float data for GLODAP only test
if glodap_only
    glodap_idx = all_data.platform < 10^6;
    vars = fieldnames(all_data);
    for v = 1:length(vars)
        all_data.(vars{v}) = all_data.(vars{v})(glodap_idx);
    end
end
clear glodap_only glodap_idx vars v

%% create directory and file names
rfr_dir = [param1 '/Models/' dir_base];
rfr_fnames = cell(num_folds,num_clusters);
for f = 1:num_folds
    for c = 1:num_clusters
        rfr_fnames(f,c) = ...
            {['RFR_' param2 '_C' num2str(c) '_F' num2str(f) '_test']};
    end
end
kfold_dir = [param1 '/KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_name = ['RFR_output_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
fig_dir = [param1 '/Figures/KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name_1 = ['k_fold_comparison_tr' num2str(numtrees) '_lf' num2str(minLeafSize) '.png'];
fig_name_2 = ['k_fold_spatial_comparison_tr' num2str(numtrees) '_lf' num2str(minLeafSize) '.png'];

%% fit and evaluate test models (RFR)
% define model parameters
NumPredictors = ceil(sqrt(length(variables)));
% set up parallel pool
tic; parpool; fprintf('Pool initiation:'); toc;
% evaluate test models for each fold
for f = 1:num_folds
    % fit test models for each cluster
    output = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    parfor c = 1:num_clusters
      % start timing fit
      tic
      if any(all_data_clusters.clusters == c) % check for data in cluster
        % fit test model for each cluster
        RFR = ...
            fit_RFR(param2,all_data,all_data_clusters.(['c' num2str(c)]),...
            train_idx.(['f' num2str(f)]),variables,numtrees,minLeafSize,...
            NumPredictors,0,thresh);
        % stop timing fit
        fprintf(['Train RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % start timing predictions
        tic
        % predict data for each cluster
        output = ...
            run_RFR(RFR,all_data,all_data_clusters.(['c' num2str(c)]),...
            test_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing predictions
        fprintf(['Run RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model for each cluster
        if ~isfolder([pwd '/' rfr_dir]); mkdir(rfr_dir);end
        parsave([rfr_dir '/' rfr_fnames{f,c}],RFR,'RFR',output,'output');
      else
        fprintf(['Train RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('\n');
        fprintf(['Run RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('\n');
        [~]=toc;
      end
    end
end
% end parallel session
delete(gcp('nocreate'))
% calculate weighted average over each cluster using probabilities
for f = 1:num_folds
    % pre-allocate probabilities
    probs_matrix = [];
    for c = 1:num_clusters
        % assemble matrix of probabilities greater than the threshold (5%)
        probs_array = all_data_clusters.(['c' num2str(c)])(test_idx.(['f' num2str(f)]));
        probs_array(probs_array < thresh) = NaN;
        probs_matrix = [probs_matrix,probs_array];
        clear probs_array
        % load output
        load([rfr_dir '/' rfr_fnames{f,c}],'output')
        rfr_output.(['f' num2str(f)])(:,c) = output;
        clear output
    end
    rfr_output.(['f' num2str(f) '_mean']) = ...
        double(sum(rfr_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
end
% clean up
clear f c
% aggregate output from all folds
rfr_output.(['k_fold_test_' param2]) = nan(size(all_data.(param2)));
for f = 1:num_folds
    rfr_output.(['k_fold_test_' param2])(test_idx.(['f' num2str(f)])) = ...
        rfr_output.(['f' num2str(f) '_mean']);
end
% compare k-fold output to data
rfr_output.k_fold_delta = rfr_output.(['k_fold_test_' param2]) - all_data.(param2);
% calculate error stats
rfr_mean_err = mean(rfr_output.k_fold_delta);
rfr_med_err = median(rfr_output.k_fold_delta);
rfr_rmse = sqrt(mean(rfr_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'rfr_output','rfr_rmse',...
    'rfr_med_err','rfr_mean_err','-v7.3');

%% plot histogram of errors
load([kfold_dir '/' kfold_name],'rfr_output','rfr_rmse');
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.(param2) rfr_output.(['k_fold_test_' param2])],...
    'Edges',{param_edges param_edges});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([param_edges(1) param_edges(end)],[param_edges(1) param_edges(end)],'k--');
set(h,'EdgeColor','none');
xlim([param_edges(1) param_edges(end)]); ylim([param_edges(1) param_edges(end)]);
xlabel(['Measured ' param2 ' (\mumol kg^{-1})']);
ylabel(['RFR ' param2 ' (\mumol kg^{-1})']);
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text((3/5)*param_edges(end),(1/10)*param_edges(end),...
    ['RMSE = ' num2str(round(rfr_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_1]);
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
rfr_output.k_fold_delta_spatial = accumarray(subs(~idx_subs,:),...
    abs(rfr_output.k_fold_delta(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure; hold on
worldmap([-90 90],[20 380]);
setm(gca,'mapprojection','robinson');
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(lat,[lon lon(end)+1],[rfr_output.k_fold_delta_spatial ...
    rfr_output.k_fold_delta_spatial(:,end)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
caxis([0 (2/50)*param_edges(end)]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Absolute \Delta' param3];
c.FontSize = 22;
c.TickLength = 0;
mlabel off; plabel off;
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name_2]);
% clean up
clear land cmap c
close

%% clean up
clear rfr_output rfr_rmse rfr_med_err rfr_mean_err probs_matrix f c
clear num_clusters f c numtrees minLeafSize NumPredictors rfr_dir rfr_fnames
clear ans all_data all_data_clusters num_folds train_idx test_idx train_sum

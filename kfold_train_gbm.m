% k_fold_fit_test_models
%
% DESCRIPTION:
% This function uses a subset of the combined dataset to train validation
% versions of machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/15/2024

function kfold_train_gbm(param,dir_base,base_grid,file_date,...
    float_file_ext,glodap_only,num_clusters,num_folds,variables,...
    numstumps,thresh)

%% process parameter name
if strcmp(param,'o2')
    param1 = 'O2';
    param2 = 'oxygen';
    param3 = '[O_{2}]';
end

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
clear glodap_idx vars v

%% create directory and file names
gbm_dir = [param1 '/Models/' dir_base];
gbm_fnames = cell(num_folds,num_clusters);
for f = 1:num_folds
    for c = 1:num_clusters
        gbm_fnames(f,c) = ...
            {['GBM_oxygen_C' num2str(c) '_F' num2str(f) '_test']};
    end
end
kfold_dir = [param1 '/KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_name = ['GBM_output_tr' num2str(numstumps)];
fig_dir = [param1 '/Figures/KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name_1 = ['k_fold_comparison_tr' num2str(numstumps) '.png'];
fig_name_2 = ['k_fold_spatial_comparison_tr' num2str(numstumps) '.png'];

%% fit and evaluate test models (GBM)
% define model parameters

% set up parallel pool
tic; parpool(6); fprintf('Pool initiation:'); toc;
% fit test models for each fold
for f = 1:num_folds
    % fit test models for each cluster
    output = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    parfor c = 1:num_clusters
      % start timing fit
      tic
      if any(all_data_clusters.clusters == c) % check for data in cluster
        % fit test model for each cluster
        GBM = ...
            fit_GBM('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
            train_idx.(['f' num2str(f)]),variables,numstumps,thresh);
        % stop timing fit
        fprintf(['Train GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % start timing predictions
        tic
        % predict data for each cluster
        output = ...
            run_GBM(GBM,all_data,all_data_clusters.(['c' num2str(c)]),...
            test_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing predictions
        fprintf(['Run GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model for each cluster
        if ~isfolder([pwd '/' gbm_dir]); mkdir(gbm_dir); end
        parsave([gbm_dir '/' gbm_fnames{f,c}],GBM,'GBM',output,'output');
      else
        fprintf(['Train GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('\n');
        fprintf(['Run GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('\n');
        [~]=toc;
      end
    end
end
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
        load([gbm_dir '/' gbm_fnames{f,c}],'output')
        gbm_output.(['f' num2str(f)])(:,c) = output;
        clear output
    end
    gbm_output.(['f' num2str(f) '_mean']) = ...
        double(sum(gbm_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
end
% clean up
clear f c
% aggregate output from all folds
gbm_output.k_fold_test_oxygen = nan(size(all_data.oxygen));
for f = 1:num_folds
    gbm_output.k_fold_test_oxygen(test_idx.(['f' num2str(f)])) = ...
        gbm_output.(['f' num2str(f) '_mean']);
end
% compare k-fold output to data
gbm_output.k_fold_delta = gbm_output.k_fold_test_oxygen - all_data.oxygen;
% calculate error stats
gbm_mean_err = mean(gbm_output.k_fold_delta);
gbm_med_err = median(gbm_output.k_fold_delta);
gbm_rmse = sqrt(mean(gbm_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'gbm_output','gbm_rmse',...
    'gbm_med_err','gbm_mean_err','-v7.3');
clear gbm_output gbm_rmse gbm_med_err gbm_mean_err

%% plot histogram of errors
load([kfold_dir '/' kfold_name],'gbm_output','gbm_rmse');
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.oxygen gbm_output.k_fold_test_oxygen],...
    'Edges',{0:5:500 0:5:500});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([0 500],[0 500],'k--');
set(h,'EdgeColor','none');
xlim([0 500]); ylim([0 500]);
xlabel('Measured Oxygen (\mumol kg^{-1})');
ylabel('Gradient-Boosted Machine Oxygen (\mumol kg^{-1})');
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(300,50,['RMSE = ' num2str(round(gbm_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
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
gbm_output.k_fold_delta_spatial = accumarray(subs(~idx_subs,:),...
    abs(gbm_output.k_fold_delta(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure; hold on
worldmap([-90 90],[20 380]);
setm(gca,'mapprojection','robinson');
set(gcf,'units','inches','position',[0 5 20 10]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(lat,[lon lon(end)+1],[gbm_output.k_fold_delta_spatial ...
    gbm_output.k_fold_delta_spatial(:,end)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
caxis([0 20]);
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
clear gbm_output gbm_rmse gbm_med_err gbm_mean_err probs_matrix
clear num_clusters numtrees minLeafSize NumPredictors gbm_dir gbm_fnames
clear all_data all_data_clusters train_idx test_idx train_sum

end
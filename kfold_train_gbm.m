% k_fold_fit_test_models
%
% DESCRIPTION:
% This function uses a subset of the combined dataset to train validation
% versions of machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% load data cluster indices
load(['Data/k_fold_data_indices_'  base_grid '_' num2str(num_clusters) ...
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
gbm_dir = ['Models/' base_grid '/GBM/GBM_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numstumps)];
gbm_fnames = cell(num_folds,num_clusters);
for f = 1:num_folds
    for c = 1:num_clusters
        gbm_fnames(f,c) = ...
            {['GBM_oxygen_C' num2str(c) '_F' num2str(f) '_test']};
    end
end
kfold_dir = ['KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_name = ['GBM_output_tr' num2str(numstumps)];
fig_dir = ['Figures/KFold/GBM/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name = ['k_fold_comparison_tr' num2str(numstumps) '.png'];

%% fit and evaluate test models (GBM)
% define model parameters

% fit test models for each fold
% LSBoost cannot run in parallel
for f = 1:num_folds
    % fit test models for each cluster
    gbm_output.(['f' num2str(f)]) = ...
        nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    for c = 1:num_clusters
      if any(all_data_clusters.clusters == c) % check for data in cluster
        % start timing fit
        tic
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
        gbm_output.(['f' num2str(f)])(:,c) = ...
            run_GBM(GBM,all_data,all_data_clusters.(['c' num2str(c)]),...
            test_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing predictions
        fprintf(['Run GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model for each cluster
        if ~isfolder([pwd '/' gbm_dir]); mkdir(gbm_dir);end
        save([gbm_dir '/' gbm_fnames{f,c}],'GBM','-v7.3');
        % clean up
        clear GBM
      else
        fprintf(['Train GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('');
        fprintf(['Run GBM - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        [~]=toc;
      end
    end
    % assemble matrix of probabilities greater than the threshold (5%)
    probs_matrix = [];
    for c = 1:num_clusters
        probs_array = all_data_clusters.(['c' num2str(c)])(test_idx.(['f' num2str(f)]));
        probs_array(probs_array < thresh) = NaN;
        probs_matrix = [probs_matrix,probs_array];
        clear probs_array
    end
    % calculate weighted average over each cluster using probabilities
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

%% plot histogram of errors
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
exportgraphics(gcf,[fig_dir '/' fig_name]);
% clean up
clear counts bin_centers h p myColorMap
close

%% clean up
clear gbm_output gbm_rmse gbm_med_err gbm_mean_err probs_matrix
clear num_clusters numtrees minLeafSize NumPredictors gbm_dir gbm_fnames
clear all_data all_data_clusters train_idx test_idx train_sum

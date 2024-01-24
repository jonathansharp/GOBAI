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
clear glodap_only glodap_idx vars v

%% create directory and file names
rfr_dir = ['Models/' base_grid '/RFR/RFR_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
rfr_fnames = cell(num_folds,num_clusters);
for f = 1:num_folds
    for c = 1:num_clusters
        rfr_fnames(f,c) = ...
            {['RFR_oxygen_C' num2str(c) '_F' num2str(f) '_test']};
    end
end
kfold_dir = ['KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_name = ['RFR_output_tr' num2str(numtrees) '_lf' num2str(minLeafSize)];
fig_dir = ['Figures/KFold/RFR/' base_grid '_c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name = ['k_fold_comparison_tr' num2str(numtrees) '_lf' num2str(minLeafSize) '.png'];

%% fit and evaluate test models (RFR)
% define model parameters
NumPredictors = ceil(sqrt(length(variables)));
% set up parallel pool
parpool;
% evaluate test models for each fold
for f = 1:num_folds
    % pre-allocate output for each fold
    rfr_output.(['f' num2str(f)]) = ...
        nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    for c = 1:num_clusters
      if any(all_data_clusters.clusters == c) % check for data in cluster
        % start timing fit
        tic
        % fit test model for each cluster
        RFR = ...
            fit_RFR('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
            train_idx.(['f' num2str(f)]),variables,numtrees,minLeafSize,...
            NumPredictors,0,thresh);
        % stop timing fit
        fprintf(['Train RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % start timing predictions
        tic
        % predict data for each cluster
        rfr_output.(['f' num2str(f)])(:,c) = ...
            run_RFR(RFR,all_data,all_data_clusters.(['c' num2str(c)]),...
            test_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing predictions
        fprintf(['Run RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model for each cluster
        if ~isfolder([pwd '/' rfr_dir]); mkdir(rfr_dir);end
        save([rfr_dir '/' rfr_fnames{f,c}],'RFR','-v7.3');
        % clean up
        clear RFR
      else
        fprintf(['Train RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
        fprintf('');
        fprintf(['Run RFR - Fold #' num2str(f) ', Cluster #' num2str(c) ': N/A']);
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
    rfr_output.(['f' num2str(f) '_mean']) = ...
        double(sum(rfr_output.(['f' num2str(f)]).*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    clear probs_matrix c
end
% end parallel session
delete(gcp('nocreate'))
% aggregate output from all folds
rfr_output.k_fold_test_oxygen = nan(size(all_data.oxygen));
for f = 1:num_folds
    rfr_output.k_fold_test_oxygen(test_idx.(['f' num2str(f)])) = ...
        rfr_output.(['f' num2str(f) '_mean']);
end
% compare k-fold output to data
rfr_output.k_fold_delta = rfr_output.k_fold_test_oxygen - all_data.oxygen;
% calculate error stats
rfr_mean_err = mean(rfr_output.k_fold_delta);
rfr_med_err = median(rfr_output.k_fold_delta);
rfr_rmse = sqrt(mean(rfr_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'rfr_output','rfr_rmse',...
    'rfr_med_err','rfr_mean_err','-v7.3');

%% plot histogram of errors
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.oxygen rfr_output.k_fold_test_oxygen],...
    'Edges',{0:5:500 0:5:500});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([0 500],[0 500],'k--');
set(h,'EdgeColor','none');
xlim([0 500]); ylim([0 500]);
xlabel('Measured Oxygen (\mumol kg^{-1})');
ylabel('RFR Oxygen (\mumol kg^{-1})');
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(300,50,['RMSE = ' num2str(round(rfr_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name]);
% clean up
clear counts bin_centers h p myColorMap
close

%% clean up
clear rfr_output rfr_rmse rfr_med_err rfr_mean_err probs_matrix f c
clear num_clusters f c numtrees minLeafSize NumPredictors rfr_dir rfr_fnames
clear ans all_data all_data_clusters num_folds train_idx test_idx train_sum

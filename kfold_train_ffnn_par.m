% k_fold_fit_test_models
%
% DESCRIPTION:
% This function uses a subset of the combined dataset to train validation
% versions of machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/17/2023

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% load data cluster indices
load(['Data/k_fold_data_indices_'  num2str(num_clusters) '_' num2str(num_folds) '_'...
    file_date float_file_ext '.mat'],'num_folds','train_idx','test_idx');

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
ffnn_dir = ['Models/FFNN/FFNN_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
ffnn_fnames = cell(5,10); % for models
ffnn_onames = cell(5,10); % for output
for f = 1:num_folds
    for c = 1:num_clusters
        ffnn_fnames(f,c) = ...
            {['FFNN_oxygen_C' num2str(c) '_F' num2str(f) '_test']};
        ffnn_onames(f,c) = ...
            {['FFNN_oxygen_C' num2str(c) '_F' num2str(f) '_test_output']};
    end
end
kfold_dir = ['KFold/FFNN/c' num2str(num_clusters) '_' file_date float_file_ext];
kfold_name = ['FFNN_output_train' num2str(100*train_ratio) '_val' num2str(100*val_ratio) ...
    '_test' num2str(100*val_ratio)];
fig_dir = ['Figures/KFold/FFNN/c' num2str(num_clusters) '_' file_date float_file_ext];
fig_name = ['k_fold_comparison_train' num2str(100*train_ratio) '_val' num2str(100*val_ratio) ...
    '_test' num2str(100*val_ratio) '.png'];

%% fit and evaluate test models (NN)
% define model parameters
nodes1 = [5 10 15];
nodes2 = [15 10 5];
% set up parallel pool
parpool;
% fit and evaluate test models for each fold
parfor fc = 1:num_folds*num_clusters
    % fit test models for each cluster
    %ffnn_output.(['f' num2str(f)]) = ...
    output = ...
        nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    parfor c = 1:num_clusters
        % start timing fit
        tic
        % fit test model for each cluster
        FFNN = ...
            fit_FFNN('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
            train_idx.(['f' num2str(f)]),variables,nodes1,nodes2,...
            train_ratio,val_ratio,test_ratio,thresh);
        % stop timing fit
        fprintf(['Train FFNN - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % start timing predictions
        tic
        % predict data for each cluster
        %ffnn_output.(['f' num2str(f)])(:,c) = ...
        output = ...
            run_FFNN(FFNN,all_data,all_data_clusters.(['c' num2str(c)]),...
            test_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing predictions
        fprintf(['Run FFNN - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model and output for each cluster
        if ~isfolder([pwd '/' ffnn_dir]); mkdir(ffnn_dir); end
        parsave_v1([ffnn_dir '/' ffnn_fnames{f,c}],FFNN);
        parsave_v1([ffnn_dir '/' ffnn_onames{f,c}],output);
        % clean up
        %clear FFNN output
    end
end
% end parallel session
delete(gcp('nocreate'));
for f = 1:num_folds
    % assemble matrix of probabilities greater than the threshold (5%)
    probs_matrix = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    output_matrix = nan(sum(test_idx.(['f' num2str(f)])),num_clusters);
    for c = 1:num_clusters
        % assemble matrix of probabilities
        probs_array = all_data_clusters.(['c' num2str(c)])(test_idx.(['f' num2str(f)]));
        probs_array(probs_array < thresh) = NaN;
        probs_matrix(:,c) = probs_array;
        clear probs_array
        % load and assemble matrix of test outputs
        output = load([ffnn_dir '/' ffnn_onames{f,c}],'output');
        output_matrix(:,c) = output;
        clear output
    end
    % calculate weighted average over each cluster using probabilities
    %ffnn_output.(['f' num2str(f) '_mean']) = ...
    output_mean = ...
        double(sum(output_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    save([ffnn_dir '/FFNN_oxygen_F' num2str(f) '_test_output_mean'],'output_mean');
    clear output_matrix probs_matrix output_mean c
end
% aggregate output from all folds
ffnn_output.k_fold_test_oxygen = nan(size(all_data.oxygen));
for f = 1:num_folds
    load([ffnn_dir '/FFNN_oxygen_F' num2str(f) '_test_output_mean'],'output_mean');
    ffnn_output.k_fold_test_oxygen(test_idx.(['f' num2str(f)])) = output_mean;
    clear output_mean
end
% compare k-fold output to data
ffnn_output.k_fold_delta = ffnn_output.k_fold_test_oxygen - all_data.oxygen;
% calculate error stats
ffnn_mean_err = mean(ffnn_output.k_fold_delta);
ffnn_med_err = median(ffnn_output.k_fold_delta);
ffnn_rmse = sqrt(mean(ffnn_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/' kfold_dir]); mkdir(kfold_dir); end
save([kfold_dir '/' kfold_name],'ffnn_output','ffnn_rmse',...
    'ffnn_med_err','ffnn_mean_err','-v7.3');
clear ffnn_output ffnn_rmse ffnn_med_err ffnn_mean_err

%% plot histogram of errors
load([kfold_dir '/' kfold_name],'ffnn_output','ffnn_rmse');
figure('visible','off'); hold on;
set(gca,'fontsize',12);
set(gcf,'position',[100 100 600 400]);
[counts,bin_centers] = hist3([all_data.oxygen ffnn_output.k_fold_test_oxygen],...
    'Edges',{0:5:500 0:5:500});
h=pcolor(bin_centers{1},bin_centers{2},counts');
plot([0 500],[0 500],'k--');
set(h,'EdgeColor','none');
xlim([0 500]); ylim([0 500]);
xlabel('Measured Oxygen (\mumol kg^{-1})');
ylabel('FFNN Oxygen (\mumol kg^{-1})');
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(300,50,['RMSE = ' num2str(round(ffnn_rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/' fig_dir]); mkdir(fig_dir); end
exportgraphics(gcf,[fig_dir '/' fig_name]);
% clean up
clear counts bin_centers h p myColorMap
close

%% clean up
clear ffnn_output ffnn_rmse ffnn_med_err ffnn_mean_err probs_matrix f c
clear f c NumPredictors ffnn_dir ffnn_fnames ffnn_onames ans all_data
clear all_data_clusters train_idx test_idx train_sum fig_dir fig_name

% k_fold_evaluate_test_models
%
% DESCRIPTION:
% This function uses validation versions of machine learning models to
% evaluate a reserved subset of data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load combined data
load_interpolated_combined_data_to_workspace

%% load data clusters
load_gmm_data_clusters

%% load k-fold evaluation indices
load_k_fold_data_indices

%% set up variables and clusters
% define variables to use
variables = define_variables_for_model();
% process cluster names
clusters = fieldnames(all_data_clusters); % obtain list of clusters
clusters(1) = []; % remove cluster identifiers from list
numClusts = length(clusters); % define number of clusters
% define cluster probability threshold
thresh = 0.05;

%% average across RFR, FFNN, SVM
% load RFR results
load(['KFold/RFR_output_'  datestr(date) '.mat'],'rfr_output');
% load FFNN results
load(['KFold/FFNN_output_'  datestr(date) '.mat'],'ffnn_output');
% load SVM results
load(['KFold/SVM_output_'  datestr(date) '.mat'],'svm_output');
% average across models
ens_output.k_fold_test_oxygen = ...
    mean([rfr_output.k_fold_test_oxygen ffnn_output.k_fold_test_oxygen ...
    svm_output.k_fold_test_oxygen],2);
% clean up
clear rfr_output ffnn_output svm_output
% compare k-fold output to data
ens_output.k_fold_delta = ens_output.k_fold_test_oxygen - all_data.oxygen;
% calculate error stats
ens_mean_err = mean(ens_output.k_fold_delta);
ens_med_err = median(ens_output.k_fold_delta);
ens_rmse = sqrt(mean(ens_output.k_fold_delta.^2));
% save predicted data
if ~isfolder([pwd '/KFold/ENS_output']); mkdir('KFold/ENS_output'); end
save(['KFold/ENS_output_'  datestr(date) '.mat'],...
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
ylabel('Ensemble Model Oxygen (\mumol kg^{-1})');
myColorMap = flipud(hot(256.*32));
myColorMap(1,:) = 1;
colormap(myColorMap);
set(gca,'ColorScale','log');
caxis([1e0 1e5]);
c=colorbar;
c.Label.String = 'log_{10}(Bin Counts)';
text(300,50,['RMSE = ' num2str(round(rmse,1)) '\mumol kg^{-1}'],'fontsize',12);
if ~isfolder([pwd '/Figures/KFold']); mkdir('KFold/KFold'); end
exportgraphics(gcf,'Figures/KFold/ENS_k_fold_comparison.png');

%% clean up
clear

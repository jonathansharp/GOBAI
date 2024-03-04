% train_rfr
%
% DESCRIPTION:
% This function uses the combined dataset to train
% random forest regressions in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/3/2024

%% load combined data
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

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
rfr_dir = ['Models/' dir_base];
rfr_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    rfr_fnames(c) = ...
        {['RFR_oxygen_C' num2str(c)]};
end
gobai_rfr_dir = ...
    ['Data/GOBAI/' base_grid '/RFR/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];

%% fit RFRs using all data

% start timing training
tic

% set up parallel pool
p = setup_pool(numWorkers_train);

% define model parameters
NumPredictors = ceil(sqrt(length(variables)));

% fit models for each cluster
parfor c = 1:num_clusters

  % start timing fit
  tic

  % check for data in cluster
  if any(all_data_clusters.clusters == c) 

    % fit model for each cluster
    RFR = ...
        fit_RFR('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
        true(size(all_data.platform)),variables,numtrees,minLeafSize,...
            NumPredictors,0,thresh);

    % save model for each cluster
    if ~isfolder([pwd '/' rfr_dir]); mkdir(rfr_dir);end
    parsave([rfr_dir '/' rfr_fnames{c}],RFR,'RFR');

    % stop timing fit
    fprintf(['Train RFR - Cluster #' num2str(c) ': ']);
    toc

  else

    % stop timing fit
    fprintf(['Train RFR - Cluster #' num2str(c) ': N/A']);
    fprintf('\n');
    [~]=toc;

  end

end

% clean up
clear all_data all_data_clusters

% end parallel session
delete(gcp('nocreate'));

% stop timing training
fprintf('RFR Training: ');
toc

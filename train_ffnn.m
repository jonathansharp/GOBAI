% train_ffnn
%
% DESCRIPTION:
% This function uses the combined dataset to train
% neural networks in GMM clusters.
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
ffnn_dir = ['Models/' dir_base];
ffnn_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    ffnn_fnames(c) = ...
        {['FFNN_oxygen_C' num2str(c)]};
end
gobai_ffnn_dir = ...
    ['Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];

%% fit FFNNs using all data

% set up parallel pool
p = setup_pool(numWorkers_train);

% start timing training
tic

% define model parameters
nodes1 = [5 10 15];
nodes2 = [15 10 5];

% fit models for each cluster
parfor c = 1:num_clusters

  % start timing fit
  tic

  % check for data in cluster
  if any(all_data_clusters.clusters == c)
    
    % fit model for each cluster
    FFNN = ...
        fit_FFNN('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
        true(size(all_data.platform)),variables,nodes1,nodes2,...
        train_ratio,val_ratio,test_ratio,thresh);

    % save model for each cluster
    if ~isfolder([pwd '/' ffnn_dir]); mkdir(ffnn_dir);end
    % save([ffnn_dir '/' ffnn_fnames{c}],'FFNN','-v7.3');
    parsave([ffnn_dir '/' ffnn_fnames{c}],FFNN,'FFNN');

    % clean up
    % clear FFNN

    % stop timing fit
    fprintf(['Train FFNN - Cluster #' num2str(c) ': ']);
    toc

  else

    % stop timing fit
    fprintf(['Train FFNN - Cluster #' num2str(c) ': N/A']);
    fprintf('\n');
    [~]=toc;

  end

end

% clean up
clear all_data all_data_clusters

% end parallel session
delete(gcp('nocreate'));

% stop timing training
fprintf('FFNN Training: ');
toc

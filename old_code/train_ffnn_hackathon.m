% train_ffnn_hackathon
%
% DESCRIPTION:
% This function uses the combined dataset to train
% neural networks in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/13/2024

%% initiate profile
profile on

%% load configuration parameters
gobai_o2_initiate;
load_standard_config_files;
load('Config/base_config_RFROM.mat'); % base grid
load('Config/predict_years_config_04.mat'); % only for 2004
dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
    float_file_ext;train_ratio;val_ratio;test_ratio}); % directory name base
fpath = '/raid'; % for RFROM
% fpath = pwd; % for RG

%% load combined data (created by gobai_o2_load.m)
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters (created by gobai_o2_cluster_data.m)
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

% start timing training
tic

% define model parameters
nodes1 = [2];
nodes2 = [2];

% set up parallel pool
tic; parpool(8); fprintf('Pool initiation:'); toc;

% fit models for just the first cluster
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

%% end and save profile
p=profile('info');
profsave(p,'profiles/train_ffnn_hackathon')
profile off

% train_rfr_hackathon
%
% DESCRIPTION:
% This function uses the combined dataset to train
% random forest regressions in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/3/2024

%% initiate profile
% profile on

%% load configuration parameters
gobai_o2_initiate;
load_standard_config_files;
numtrees = 5;
load('Config/base_config_RFROM.mat'); % base grid
load('Config/predict_years_config_04.mat'); % only for 2004
dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
fpath = '/raid'; % for RFROM
% fpath = pwd; % for RG

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

% define model parameters
NumPredictors = ceil(sqrt(length(variables)));

% set up parallel pool
%tic; parpool(8); fprintf('Pool initiation:'); toc;

% fit models for just the first cluster
c = 1;

  % check for data in cluster
  if any(all_data_clusters.clusters == c) 

    % start timing fit
    tic

    % normalize data
    

    % fit model for each cluster
    RFR = ...
        fit_RFR('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
        true(size(all_data.platform)),variables,numtrees,minLeafSize,...
            NumPredictors,0,thresh);

    % save model for each cluster
    if ~isfolder([pwd '/' rfr_dir]); mkdir(rfr_dir);end
    save([rfr_dir '/' rfr_fnames{c}],'RFR','-v7.3');

    % clean up
    clear RFR

    % stop timing fit
    fprintf(['Train RFR - Cluster #' num2str(c) ': ']);
    toc

  else

    % stop timing fit
    fprintf(['Train RFR - Cluster #' num2str(c) ': N/A']);
    fprintf('\n');
    [~]=toc;

  end


% clean up
clear all_data all_data_clusters

% end parallel session
delete(gcp('nocreate'));

% stop timing training
fprintf('RFR Training: ');
toc

%% end and save profile
% p=profile('info');
% profsave(p,'profiles/train_rfr_hackathon')
% profile off

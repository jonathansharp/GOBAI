% k_fold_split_data
%
% DESCRIPTION:
% This function splits the data into k-fold subsets.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% remove float data for GLODAP only test
if glodap_only
    glodap_idx = all_data.platform < 10^6;
    vars = fieldnames(all_data);
    for v = 1:length(vars)
        all_data.(vars{v}) = all_data.(vars{v})(glodap_idx);
    end
end
clear glodap_only glodap_idx vars v

%% split data into testing and training sets equal to k
rng(8); % for reproducibility
plat_ids = unique(all_data.platform); % unique platform IDs
num_plats = length(plat_ids); % number of unique platforms
ints = randperm(num_plats)'; % assign sequential integers to platforms
for f = 1:numFolds
    % determine index of platforms that belong to fold
    split_idx = ...
        ints > (f-1) * (num_plats/numFolds) & ints <= f * (num_plats/numFolds);
    % determine platforms that belong to fold
    test_plat.(['f' num2str(f)]) = plat_ids(split_idx);
    train_plat.(['f' num2str(f)]) = plat_ids(~split_idx);
    % determine index of data from platforms that belong to fold
    test_idx.(['f' num2str(f)]) = ...
        ismember(all_data.platform,test_plat.(['f' num2str(f)]));
    train_idx.(['f' num2str(f)]) = ...
        ismember(all_data.platform,train_plat.(['f' num2str(f)]));
end
clear plat_ids num_plats ints split_idx test_plat train_plat

%% process cluster names
clusters = fieldnames(all_data_clusters); % obtain list of clusters
clusters(1) = []; % remove cluster identifiers from list
numClusts = length(clusters); % define number of clusters

%% determine and save number of training points for each cluster and fold
% pre-allocate matrix for number of test/training data points
training_data_points = nan(numFolds,numClusts);
test_data_points = nan(numFolds,numClusts);
% determine number of test/training points for each cluster and fold
for f = 1:numFolds
    for c = 1:numClusts
        idx_train = train_idx.(['f' num2str(f)]) & ...
            all_data_clusters.(['c' num2str(c)]) > thresh;
        training_data_points(f,c) = sum(idx_train(:));
        idx_test = test_idx.(['f' num2str(f)]) & ...
            all_data_clusters.(['c' num2str(c)]) > thresh;
        test_data_points(f,c) = sum(idx_test(:));
    end
end
% save tables of test/training data points
training_data_points_table = array2table(training_data_points,...
    'RowNames',sprintfc('%d',1:numFolds),'VariableNames',clusters);
test_data_points_table = array2table(test_data_points,...
    'RowNames',sprintfc('%d',1:numFolds),'VariableNames',clusters);
if ~isfolder([pwd '/Data']); mkdir('Data'); end
writetable(training_data_points_table,['Data/training_data_points_' datestr(date) '.csv']);
writetable(test_data_points_table,['Data/test_data_points_' datestr(date) '.csv']);
clear training_data_points training_data_points_table idx_train
clear test_data_points test_data_points_table idx_test

%% save k-fold evaluation indices
if ~isfolder([pwd '/Data']); mkdir('Data'); end
save(['Data/k_fold_data_indices_'  num2str(num_clusters) '_' num2str(numFolds) '_'...
    file_date float_file_ext '.mat'],'numFolds','train_idx','test_idx','-v7.3');

% k_fold_split_data
%
% DESCRIPTION:
% This function splits the data into k-fold subsets.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/18/2023

function kfold_split_data(param,base_grid,file_date,float_file_ext,...
    glodap_only,num_clusters,num_folds,thresh)

%% process parameter name
param1 = param_name(param);

if ~exist([param1 '/Data/k_fold_data_indices_'  base_grid '_' num2str(num_clusters) ...
        '_' num2str(num_folds) '_' file_date float_file_ext '.mat'],'file')

%% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load([param1 '/Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
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
for f = 1:num_folds
    % determine index of platforms that belong to fold
    split_idx = ...
        ints > (f-1) * (num_plats/num_folds) & ints <= f * (num_plats/num_folds);
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
training_data_points = nan(num_folds,numClusts);
test_data_points = nan(num_folds,numClusts);
% determine number of test/training points for each cluster and fold
for f = 1:num_folds
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
    'RowNames',sprintfc('%d',1:num_folds),'VariableNames',clusters);
test_data_points_table = array2table(test_data_points,...
    'RowNames',sprintfc('%d',1:num_folds),'VariableNames',clusters);
if ~isfolder([pwd '/' param1 '/Data']); mkdir([param1 '/Data']); end
writetable(training_data_points_table,[param1 '/Data/training_data_points_' datestr(date) '.csv']);
writetable(test_data_points_table,[param1 '/Data/test_data_points_' datestr(date) '.csv']);
clear training_data_points training_data_points_table idx_train
clear test_data_points test_data_points_table idx_test

%% save k-fold evaluation indices
if ~isfolder([pwd '/' param1 '/Data']); mkdir([param1 '/Data']); end
save([param1 '/Data/k_fold_data_indices_'  base_grid '_' num2str(num_clusters) '_' num2str(num_folds) '_'...
    file_date float_file_ext '.mat'],'num_folds','train_idx','test_idx','-v7.3');

end

end
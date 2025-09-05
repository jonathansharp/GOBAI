% k_fold_split_data
%
% DESCRIPTION:
% This function splits the data into k-fold subsets.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/20/2025

function kfold_split_data(param_props,file_date,float_file_ext,...
    num_clusters,num_folds,thresh,flt,gld,ctd)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% load combined data
load([param_props.dir_name '/Data/processed_all_' param_props.file_name ...
    '_data_' float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load([param_props.dir_name '/Data/GMM_' num2str(num_clusters) '/all_data_clusters_' num2str(num_clusters) '_' ...
    float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],'all_data_clusters');

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
if ~isfolder([param_props.dir_name '/Data']); mkdir([param_props.dir_name '/Data']); end
writetable(training_data_points_table,[param_props.dir_name '/Data/training_data_points_' datestr(date) '.csv']);
writetable(test_data_points_table,[param_props.dir_name '/Data/test_data_points_' datestr(date) '.csv']);
clear training_data_points training_data_points_table idx_train
clear test_data_points test_data_points_table idx_test

%% save k-fold evaluation indices
if ~isfolder([param_props.dir_name '/Data']); mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/k_fold_data_indices_' num2str(num_clusters) ...
    '_' num2str(num_folds) '_' float_ext glodap_ext ctd_ext '_' ...
    file_date float_file_ext '.mat'],'num_folds','train_idx','test_idx','-v7.3');

%% display information
disp(['Data split into ' num2str(num_folds) ' folds']);

end
% k_fold_fit_test_models
%
% DESCRIPTION:
% This function uses a subset of the combined dataset to train validation
% versions of machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% **EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
glodap_only = false;

%% load combined data
load_interpolated_combined_data_to_workspace

%% load data clusters
load_gmm_data_clusters

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
numFolds = 5; % number of folds
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
% save k-fold evaluation indices
if ~isfolder([pwd '/Data']); mkdir('Data'); end
save(['Data/k_fold_data_indices_' num2str(numFolds) '.mat'],...
    'numFolds','train_idx','test_idx','-v7.3');

%% set up variables and clusters
% define variables to use
variables = define_variables_for_model();
% process cluster names
clusters = fieldnames(all_data_clusters); % obtain list of clusters
clusters(1) = []; % remove cluster identifiers from list
numClusts = length(clusters); % define number of clusters
% define cluster probability threshold
thresh = 0.05;

%% determine and save number of training/test points for each cluster and fold
% pre-allocate matrix for number of training data points
training_data_points = nan(numFolds,numClusts);
% determine number of training points for each cluster and fold
for f = 1:numFolds
    for c = 1:numClusts
        idx_train = train_idx.(['f' num2str(f)]) & ...
            all_data_clusters.(['c' num2str(c)]) > thresh;
        training_data_points(f,c) = sum(idx_train(:));
    end
end
% save table of training data points
training_data_points_table = array2table(training_data_points,...
    'RowNames',sprintfc('%d',1:numFolds),'VariableNames',clusters);
if ~isfolder([pwd '/Data']); mkdir('Data'); end
writetable(training_data_points_table,['Data/training_data_points_' datestr(date) '.csv']);
clear training_data_points training_data_points_table idx_train

%% fit test models (RFR)
% set up parallel pool if none exists
p = gcp('nocreate'); if isempty(p); parpool; end
% fit test models for each fold
for f = 1:numFolds
    % fit test models for each cluster
    for c = 1:numClusts
        % define model parameters
        numtrees = 100;
        minLeafSize = 5;
        NumPredictors = ceil(sqrt(length(variables)));
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
        % save test model for each cluster
        if ~isfolder([pwd '/Models/RFR']); mkdir('Models/RFR'); end
        save(['Models/RFR/RFR_oxygen_C' num2str(c) '_F' num2str(f) '_test'],'RFR','-v7.3');
        % clean up
        clear RFR
    end
end
% end parallel session
delete(gcp('nocreate'))

%% fit test models (FFNN)
% set up parallel pool if none exists
p = gcp('nocreate'); if isempty(p); parpool; end
% fit test models for each fold
for f = 1:numFolds
    % fit test models for each cluster
    for c = 1:numClusts
        % define model parameters
        nodes1 = [10 15 20];
        nodes2 = [20 15 10];
        train_ratio = 0.7;
        val_ratio = 0.15;
        test_ratio = 0.15;
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
        % save test model for each cluster
        if ~isfolder([pwd '/Models/FFNN']); mkdir('Models/FFNN'); end
        save(['Models/FFNN/FFNN_oxygen_C' num2str(c) '_F' num2str(f) '_test'],'FFNN','-v7.3');
        % clean up
        clear FFNN
    end
end
% end parallel session
delete(gcp('nocreate'))

%% fit test models (SVM)
% fit test models for each fold
for f = 1:numFolds
    % fit test models for each cluster
    for c = 1:numClusts
        % define model parameters

        % start timing fit
        tic
        % fit test model for each cluster
        SVM = ...
            fit_SVM('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
            train_idx.(['f' num2str(f)]),variables,thresh);
        % stop timing fit
        fprintf(['Train SVM - Fold #' num2str(f) ', Cluster #' num2str(c) ': ']);
        toc
        % save test model for each cluster
        if ~isfolder([pwd '/Models/SVM']); mkdir('Models/SVM'); end
        save(['Models/SVM/SVM_oxygen_C' num2str(c) '_F' num2str(f) '_test'],'SVM','-v7.3');
        % clean up
        clear SVM
    end
end

%% clean up
clear variables clusters numClusts f c numtrees minLeafSize NumPredictors
clear ans all_data all_data_clusters numFolds train_idx test_idx train_sum

% create directory
if ~isfolder('Config'); mkdir('Config'); end

%% create data download configuration file for adjusted and DMQC
clear
snap_date = 202309;
glodap_year = 2022;
data_modes = {'A' 'D'};
float_file_ext = [];
for m = 1:length(data_modes)
    float_file_ext = [float_file_ext '_' data_modes{m}];
end
snap_download = 0; 
save(['Config/load_data_config' float_file_ext '.mat']);
clear

%% create data download configuration file for DMQC only
clear
snap_date = 202309;
glodap_year = 2022;
data_modes = {'D'};
float_file_ext = [];
for m = 1:length(data_modes)
    float_file_ext = [float_file_ext '_' data_modes{m}];
end
snap_download = 0; 
save(['Config/load_data_config' float_file_ext '.mat']);
clear


%% create cluster configuration files
clear
for num_clusters = 2:1:30
    save(['Config/cluster_config_' num2str(num_clusters) '.mat']);
end
clear

%% create k-fold configuration files
clear
glodap_only = false; % EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
thresh = 0.05;
for numFolds = [5 10] % number of folds
    for numClusts = [5 10 15 20 30] % number of clusters
        variables = ... % variables for algorithms
            {'latitude' 'longitude' 'pressure' 'sigma' 'temperature_cns' ...
            'salinity_abs'};
        save(['Config/kFold_config_' num2str(numClusts) '_' ...
            num2str(numFolds) '.mat']);
    end
end
clear

%% create RFR configuration files
clear
for numtrees = 100:100:1000
    for minLeafSize = [2 5 10 15 20 25 30]
        save(['Config/rfr_config_' num2str(numtrees) ...
            '_' num2str(minLeafSize) '.mat']);
    end
end
clear

%% create FFNN configuration files
clear
train_ratio = 0.7;
val_ratio = 0.15;
test_ratio = 0.15;
save(['Config/ffnn_config_' num2str(train_ratio*100) '_' ...
    num2str(val_ratio*100) '_' num2str(test_ratio*100) '.mat']);
clear
train_ratio = 0.85;
val_ratio = 0.15;
test_ratio = 0.0;
save(['Config/ffnn_config_' num2str(train_ratio*100) '_' ...
    num2str(val_ratio*100) '_' num2str(test_ratio*100) '.mat']);
clear
train_ratio = 0.8;
val_ratio = 0.1;
test_ratio = 0.1;
save(['Config/ffnn_config_' num2str(train_ratio*100) '_' ...
    num2str(val_ratio*100) '_' num2str(test_ratio*100) '.mat']);
clear

%% create GBM configuration files
clear
for numstumps = [100 200 400 600 800 1000]
    save(['Config/gbm_config_' num2str(numstumps) '.mat']);
end
clear

%% create base grid configuration files
clear
base_grid = 'RG';
save('Config/base_config_RG.mat');
clear
base_grid = 'RFROM';
save('Config/base_config_RFROM.mat');
clear

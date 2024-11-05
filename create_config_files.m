% create directory
if ~isfolder('Config'); mkdir('Config'); end

%% create systemsconfiguration files
numWorkers_train = 40;
numWorkers_predict = 20;
model_path = 'TBD';
snap_download = 0;
save('Config/hercules.mat');
clear
numWorkers_train = 2;
numWorkers_predict = 1;
model_path = '/raid/Model/CMIP6/';
snap_download = 1;
save('Config/chinook.mat');
clear

%% create data download configuration file for adjusted and DMQC
clear
snap_date = 202410;
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
data_modes = {'A' 'D'};
float_file_ext = '_A_D';
save(['Config/load_data_config' float_file_ext '.mat']);
clear

%% create data download configuration file for DMQC only
clear
snap_date = 202410;
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
glodap_year = 2023;
snap_download = 0; 
data_modes = {'D'};
float_file_ext = '_D';
save(['Config/load_data_config' float_file_ext '.mat']);
clear

%% create cluster configuration files
clear
clust_vars = {'temperature_cns' 'salinity_abs' 'pressure'};
% clust_vars = {'temperature' 'salinity' 'pressure'};
for num_clusters = 2:1:30
    save(['Config/cluster_config_' num2str(num_clusters) '.mat']);
end
clear

%% create k-fold configuration files
clear
glodap_only = false; % EDIT THIS TO 'true' TO TEST WITH GLODAP DATA ONLY
thresh = 0.05;
for num_folds = [5 10] % number of folds
    save(['Config/kfold_config_' num2str(num_folds) '.mat']);
end
clear

%% create algorithm training files

% all variables
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs' 'day_sin' 'day_cos' 'year'};
save('Config/vars_config_all.mat');

% no spatiotemporal coordinates
variables = ... % variables for algorithms
    {'sigma' 'temperature_cns' 'salinity_abs'};
save('Config/vars_config_no_coords.mat');

% no time
variables = ... % variables for algorithms
    {'latitude' 'lon_cos_1' 'lon_cos_2' 'pressure' 'sigma' ...
    'temperature_cns' 'salinity_abs'};
save('Config/vars_config_no_time.mat');

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

%% create GBM stumps configuration files
clear
for numstumps = [100 200 400 600 800 1000 2000 5000]
    save(['Config/gbm_config_stumps_' num2str(numstumps) '.mat']);
end
clear

%% create GBM bins configuration files
clear
for numbins = [10:10:100 200 500 1000]
    save(['Config/gbm_config_bins_' num2str(numbins) '.mat']);
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

%% create years to predict configuration files
% 4-year chunks
clear
y1 = 2004;
while y1+3 <= 2023
    years_to_predict = y1:y1+3;
    save(['Config/predict_years_config_' sprintf('%02d',(y1-2000)) ...
        '_' sprintf('%02d',(y1+3-2000)) '.mat'],'years_to_predict');
    y1=y1+4;
end
clear
% 1-year chunks
clear
y1 = 2004;
while y1 <= 2023
    years_to_predict = y1;
    save(['Config/predict_years_config_' sprintf('%02d',(y1-2000)) ...
        '.mat'],'years_to_predict');
    y1=y1+1;
end
clear

% assign_data_to_clusters
%
% DESCRIPTION:
% This function loads processed/combined float and glodap data and assigns
% them to clusters formed from gridded tempearture and salinity data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/31/2025

function assign_data_to_clusters(param_props,base_grid,snap_date,float_file_ext,...
    clust_vars,num_clusters,flt,gld,ctd)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% process date
date_str = num2str(snap_date);
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');

%% load combined data
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param_props.dir_name '/Data/processed_all_' ...
        param_props.file_name '_data_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.mat'],...
         'all_data','file_date');
else
    load([param_props.dir_name '/Data/' base_grid '_' ...
        param_props.file_name '_data_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.mat'],...
         'all_data','file_date');
end

%% define GMM model name
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
        num2str(num_clusters)];
    if ~isfolder(gmm_folder_name); mkdir(gmm_folder_name); end
    gmm_model_name = [gmm_folder_name '/model_' float_ext ...
        glodap_ext ctd_ext '_' date_str];
else
    gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
        base_grid '_' num2str(num_clusters)];
    if ~isfolder(gmm_folder_name); mkdir(gmm_folder_name); end
    gmm_model_name = [gmm_folder_name '/model_' float_ext ...
        glodap_ext ctd_ext '_' date_str];
end

%% assign data points and probabilities to clusters
% load GMM model
load(gmm_model_name,'gmm','C','S');
% transform to normalized arrays
predictor_matrix = [];
for v = 1:length(clust_vars)
    predictor_matrix = [predictor_matrix all_data.(clust_vars{v})];
end
X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
% assign to clusters and obtain probabilities
[all_data_clusters.clusters,~,p] = cluster(gmm,X_norm);
all_data_clusters.clusters = uint8(all_data_clusters.clusters);
% assign probabilities to data cluster structure
for c = 1:size(p,2)
    all_data_clusters.(['c' num2str(c)]) = p(:,c);
end
% save data clusters
if ~isfolder([pwd '/Data']); mkdir('Data'); end
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    save([param_props.dir_name '/Data/GMM_' num2str(num_clusters) ...
        '/all_data_clusters_' num2str(num_clusters) '_' ...
        float_ext glodap_ext ctd_ext '_' file_date float_file_ext ...
        '.mat'],'all_data_clusters','-v7.3');
else
    save([param_props.dir_name '/Data/GMM_' base_grid '_' ...
        num2str(num_clusters) '/all_data_clusters_' ...
        num2str(num_clusters) '_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.mat'],'all_data_clusters','-v7.3');
end
% display information
disp(['data assigned to ' num2str(num_clusters) ' clusters']);

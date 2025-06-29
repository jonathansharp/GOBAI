% assign_data_to_clusters
%
% DESCRIPTION:
% This function loads processed/combined float and glodap data and assigns
% them to clusters formed from gridded tempearture and salinity data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 10/05/2023

function assign_data_to_clusters(param_props,base_grid,snap_date,float_file_ext,...
    clust_vars,num_clusters)

%% process date
date_str = num2str(snap_date);
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');

%% load combined data
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    load([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
         'all_data','file_date');
else
    load([param_props.dir_name '/Data/' base_grid '_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
         'all_data','file_date');
end

%% assign data points and probabilities to clusters
% load GMM model
load([param_props.dir_name '/Data/GMM_' num2str(num_clusters) '/model_' date_str]);
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
save([param_props.dir_name '/Data/all_data_clusters_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters','-v7.3');
% display information
disp(['data assigned to ' num2str(num_clusters) ' clusters']);

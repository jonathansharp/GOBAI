% assign_data_to_clusters
%
% DESCRIPTION:
% This function loads processed/combined float and glodap data and assigns
% them to clusters formed from gridded tempearture and salinity data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 10/05/2023

function assign_data_to_clusters(param,base_grid,file_date,float_file_ext,...
    clust_vars,num_clusters)

%% process parameter name
param1 = param_name(param);

%% load combined data
load([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% assign data points and probabilities to clusters
% load GMM model
load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model']);
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
save([param1 '/Data/all_data_clusters_'  base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters','-v7.3');
clear lon_cos X_norm c C S gmm p

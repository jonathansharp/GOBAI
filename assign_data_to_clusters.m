% assign_data_to_clusters
%
% DESCRIPTION:
% This function loads processed/combined float and glodap data and assigns
% them to clusters formed from gridded tempearture and salinity data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 10/05/2023

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
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
save(['Data/all_data_clusters_'  base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters','-v7.3');
clear lon_cos X_norm c C S gmm p

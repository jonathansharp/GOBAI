% Runs all scripts to cluster data for GOBAI-O2

addpath(genpath(pwd));
try
    gmm_clustering(base_grid,clust_vars,num_clusters,numWorkers_predict);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
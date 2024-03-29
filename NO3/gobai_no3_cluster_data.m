% Runs all scripts to cluster data for GOBAI-NO3

addpath(genpath(pwd));
try
    assign_data_to_clusters('no3',base_grid,file_date,float_file_ext,...
        clust_vars,num_clusters);
    plot_data_by_cluster('no3',base_grid,file_date,float_file_ext,...
        clust_vars,num_clusters);
    kfold_split_data('no3',base_grid,file_date,float_file_ext,...
        glodap_only,num_clusters,num_folds,thresh);
    disp('success!');
    disp(['data assigned to ' num2str(num_clusters) ' clusters on ' ...
        base_grid ' grid and split into ' num2str(num_folds) ' folds']);
catch ME
    display_error_info(ME);
end
exit
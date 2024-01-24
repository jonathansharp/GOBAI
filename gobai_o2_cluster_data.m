% Runs all scripts to cluster data for GOBAI-O2

addpath(genpath(pwd));
try
    assign_data_to_clusters;
    kfold_split_data;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
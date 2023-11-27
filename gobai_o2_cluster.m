% Runs all scripts to cluster data for GOBAI-O2

addpath(genpath(pwd));
try
    gmm_clustering_rg;
    assign_data_to_clusters;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
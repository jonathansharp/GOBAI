% Runs all scripts to cluster data for GOBAI-O2

addpath(genpath(pwd));
try
    gmm_clustering;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
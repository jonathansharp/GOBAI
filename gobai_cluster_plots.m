% Runs all scripts to create cluster plots for GOBAI-O2

addpath(genpath(pwd));
try
    % plot_clusters;
    plot_cluster_animation(param,base_grid,num_clusters,pressure);
    plot_probability_animation(param,base_grid,num_clusters,pressure);
    disp('success!');
catch ME % if you can
    display_error_info(ME);
end
exit
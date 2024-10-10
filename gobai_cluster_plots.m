% Runs all scripts to create cluster plots for GOBAI-O2

try
    plot_cluster_animation(base_grid,num_clusters);
    plot_probability_animation(base_grid,num_clusters);
    disp('success!');
catch ME % if you can
    display_error_info(ME);
end
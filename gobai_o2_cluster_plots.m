% Runs all scripts to create cluster plots for GOBAI-O2

addpath(genpath(pwd));
try
    plot_data_by_cluster;
    plot_cluster_animation
    %plot_probability_animation;
    disp('success!');
catch ME % if you can
    display_error_info(ME);
end
exit
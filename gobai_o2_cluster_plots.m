% Runs all scripts to create cluster plots for GOBAI-O2

addpath(genpath(pwd));
try
    plot_cluster_animation
    plot_probability_animation;
    plot_data_by_cluster;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
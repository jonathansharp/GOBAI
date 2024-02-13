% Runs all scripts to plot GOBAI-O2 GBM gridded fields

addpath(genpath(pwd));
try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    mod_type = 'GBM';
    %plot_gobai_mean;
    plot_gobai_animation;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
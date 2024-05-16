% Runs all scripts to plot GOBAI-O2 FFNN gridded fields

addpath(genpath(pwd));
try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    mod_type = 'FFNN';
    %plot_gobai_mean;
    plot_gobai_animation('o2',dir_base,base_grid,num_clusters,mod_type)
    disp('success!');
catch ME
    display_error_info(ME);
end
%exit
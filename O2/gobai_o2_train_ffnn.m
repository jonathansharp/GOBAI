% Runs all scripts to optimize GOBAI-O2 machine learning algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    opt_rfr(base_grid,file_date,float_file_ext);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
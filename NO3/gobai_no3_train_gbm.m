% Runs all scripts to train GOBAI-NO3 GBM algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    train_gbm;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
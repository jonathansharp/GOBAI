% Runs all scripts to train GOBAI-O2 FFNN test algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    kfold_train_ffnn;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
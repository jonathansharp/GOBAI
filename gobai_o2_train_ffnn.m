% Runs all scripts to train and apply GOBAI-O2 FFNN algorithms
mpiprofile on
addpath(genpath(pwd));
try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    train_ffnn;
    disp('success!');
catch ME
    display_error_info(ME);
end
mpiprofile viewer
% exit
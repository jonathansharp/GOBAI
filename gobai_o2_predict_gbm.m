% Runs all scripts to apply GOBAI-O2 GBM algorithms to gridded fields

addpath(genpath(pwd));
try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    predict_gbm;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
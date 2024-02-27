% Runs all scripts to apply GOBAI-O2 FFNN algorithms to gridded fields

try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    predict_ffnn;
    aggregate_GOBAI_files;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
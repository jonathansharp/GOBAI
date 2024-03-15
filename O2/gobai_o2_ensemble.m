% Runs all scripts to combine GOBAI-O2 FFNN gridded fields

try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    combine_gobai;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
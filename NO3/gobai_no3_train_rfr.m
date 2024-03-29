% Runs all scripts to train and apply GOBAI-NO3 RFR algorithms

try
    numtrees=5;
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    train_rfr;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
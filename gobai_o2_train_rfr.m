% Runs all scripts to train and apply GOBAI-O2 RFR algorithms

try
    numtrees = 3;
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    train_rfr;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
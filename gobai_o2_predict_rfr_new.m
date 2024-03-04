% Runs all scripts to apply GOBAI-O2 RFR algorithms to gridded fields

try
    numtrees = 5;
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    predict_rfr_edited;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
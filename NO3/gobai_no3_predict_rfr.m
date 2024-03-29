% Runs all scripts to apply GOBAI-NO3 RFR algorithms to gridded fields

try
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    predict_rfr;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
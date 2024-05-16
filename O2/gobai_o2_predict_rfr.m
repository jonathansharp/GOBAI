% Runs all scripts to apply GOBAI-O2 RFR algorithms to gridded fields

try
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    predict_rfr('o2',dir_base,base_grid,file_date,...
        float_file_ext,num_clusters,variables,numtrees,minLeafSize,...
        thresh,numWorkers_predict,years_to_predict);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
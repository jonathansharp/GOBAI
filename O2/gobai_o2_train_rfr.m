% Runs all scripts to train and apply GOBAI-O2 RFR algorithms

try
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    train_rfr('o2',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,variables,...
        numtrees,minLeafSize,thresh);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
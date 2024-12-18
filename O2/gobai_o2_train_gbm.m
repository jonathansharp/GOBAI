% Runs all scripts to train GOBAI-O2 GBM algorithms

try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    train_gbm('o2',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,variables,...
        numstumps,numbins,thresh);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
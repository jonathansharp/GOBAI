% Runs all scripts to apply GOBAI-O2 GBM algorithms to gridded fields

try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    predict_gbm('o2',dir_base,base_grid,file_date,...
        float_file_ext,num_clusters,variables,numstumps,numbins,thresh,...
        numWorkers_predict,years_to_predict);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
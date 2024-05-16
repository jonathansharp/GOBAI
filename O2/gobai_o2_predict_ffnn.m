% Runs all scripts to apply GOBAI-O2 FFNN algorithms to gridded fields

try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    predict_ffnn('o2',dir_base,base_grid,file_date,...
        float_file_ext,num_clusters,variables,train_ratio,...
        val_ratio,test_ratio,thresh,numWorkers_predict,years_to_predict);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
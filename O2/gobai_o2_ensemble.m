% Runs all scripts to combine GOBAI-O2 FFNN gridded fields

try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    combine_gobai(param,dir_base,base_grid,file_date,float_file_ext,...
        num_clusters,variables,train_ratio,val_ratio,test_ratio,thresh,...
        numWorkers_predict,years_to_predict);
    mod_type = 'ENS';
    %plot_gobai_mean;
    plot_gobai_animation('o2',dir_base,base_grid,num_clusters,mod_type);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
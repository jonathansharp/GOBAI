% Runs all scripts to train GOBAI-NO3 FFNN test algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
    kfold_train_ffnn('no3',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,num_folds,variables,...
        train_ratio,val_ratio,test_ratio,thresh);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
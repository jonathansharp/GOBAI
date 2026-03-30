% Runs all scripts to train GOBAI-O2 GBM test algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    kfold_train_gbm('o2',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,num_folds,variables,...
        numstumps,numbins,thresh,0);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
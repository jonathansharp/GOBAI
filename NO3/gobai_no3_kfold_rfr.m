% Runs all scripts to train GOBAI-NO3 RFR test algorithms

addpath(genpath(pwd));
try
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    kfold_train_rfr('no3',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,num_folds,variables,...
        numtrees,minLeafSize,thresh);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
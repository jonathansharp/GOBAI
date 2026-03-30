% Runs all scripts to train GOBAI-O2 GBM test algorithms

load_standard_config_files;
load('Config/workers_chinook.mat');
load('Config/base_config_RFROM.mat');
numbins = [10 20 30 40 50 60 70 80 90 100];

for bv = 1:length(numbins)

try
    dir_base = create_dir_base('GBM',{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps});
    kfold_train_gbm('o2',dir_base,base_grid,file_date,...
        float_file_ext,glodap_only,num_clusters,num_folds,variables,...
        numstumps,numbins(bv),thresh,1);
    disp('success!');
catch ME
    display_error_info(ME);
end
%exit

end
% Runs all scripts to train GOBAI-O2 GBM test algorithms

addpath(genpath(pwd));
try
    kfold_train_gbm;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
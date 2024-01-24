% Runs all scripts to train GOBAI-O2 FFNN test algorithms

addpath(genpath(pwd));
try
    kfold_train_ffnn;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
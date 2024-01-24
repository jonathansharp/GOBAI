% Runs all scripts to train GOBAI-O2 RFR test algorithms

addpath(genpath(pwd));
try
    kfold_train_rfr;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
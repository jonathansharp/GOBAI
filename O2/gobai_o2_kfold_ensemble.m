% Runs all scripts to evaluate GOBAI-O2 test algorithms

addpath(genpath(pwd));
try
    kfold_ensemble;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
% Runs all scripts to evaluate GOBAI-NO3 test algorithms

addpath(genpath(pwd));
try
    kfold_ensemble;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
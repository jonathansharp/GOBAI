% Runs all scripts to split GOBAI-O2 data for testing

addpath(genpath(pwd));
try
    kfold_split_data
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
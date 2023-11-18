% Runs all scripts to split GOBAI-O2 data for testing

addpath(genpath(pwd));
try
    kfold_split_data
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
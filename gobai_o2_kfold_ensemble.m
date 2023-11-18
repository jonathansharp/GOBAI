% Runs all scripts to evaluate GOBAI-O2 test algorithms

addpath(genpath(pwd));
try
    kfold_ensemble;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
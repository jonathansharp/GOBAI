% Runs all scripts to test the GOBAI-O2 models using model output

addpath(genpath(pwd));
try
    load_GFDL_ESM4;
    % train_osse_models;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
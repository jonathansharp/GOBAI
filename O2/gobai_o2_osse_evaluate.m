% Runs all scripts to test the GOBAI-O2 models using model output

addpath(genpath(pwd));
try
    load_GFDL_ESM4(model_path);
    % train_osse_models;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
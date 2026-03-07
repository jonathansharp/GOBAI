% Runs all scripts to test the GOBAI-NO3 models using model output

addpath(genpath(pwd));
try
    load_GFDL_ESM4;
    % train_osse_models;
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
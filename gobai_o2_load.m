% Runs all scripts to import and process data for GOBAI-O2

addpath(genpath(pwd));
try
    acquire_o2_snapshot_data;
    acquire_o2_glodap_data;
    display_o2_data;
    %calculate_o2_gridding_uncertainty;
    adjust_o2_float_data;
    combine_o2_data;
    disp('success!');
catch ME
    disp(ME);
    disp(['Function: ' ME.stack(1).name]);
    disp(['Line: ' num2str(ME.stack(1).line)]);
end
exit
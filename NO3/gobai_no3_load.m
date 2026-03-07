% Runs all scripts to import and process data for GOBAI-NO3

addpath(genpath(pwd));
try
    acquire_snapshot_data('no3',data_modes,float_file_ext,snap_date,snap_download);
    acquire_glodap_data('no3',glodap_year);
    display_data('no3',float_file_ext,file_date,glodap_year);
    % calculate_gridding_uncertainty;
    adjust_no3_float_data(float_file_ext,file_date,glodap_year);
    combine_data('no3',float_file_ext,file_date,glodap_year);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
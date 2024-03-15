% Runs all scripts to import and process data for GOBAI-O2

addpath(genpath(pwd));
try
    acquire_o2_snapshot_data(data_modes,float_file_ext,snap_date,snap_download);
    acquire_o2_glodap_data(glodap_year);
    display_o2_data(float_file_ext,snap_date,glodap_year);
    %calculate_o2_gridding_uncertainty;
    adjust_o2_float_data(float_file_ext,snap_date,glodap_year);
    combine_o2_data(float_file_ext,snap_date,glodap_year);
    disp('success!');
catch ME
    display_error_info(ME);
end
exit
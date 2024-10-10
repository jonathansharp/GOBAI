% Runs all scripts to import and process data for GOBAI-O2

try
    acquire_snapshot_data('o2',data_modes,float_file_ext,snap_date,snap_download);
    acquire_glodap_data('o2',glodap_year);
    display_data('o2',float_file_ext,file_date,glodap_year);
    % calculate_gridding_uncertainty;
    adjust_o2_float_data(float_file_ext,file_date,glodap_year);
    combine_data('o2',float_file_ext,file_date,glodap_year);
    disp('success!');
catch ME
    display_error_info(ME);
end
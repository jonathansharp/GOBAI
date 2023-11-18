% load_interpolated_adjusted_data_to_workspace
%
% DESCRIPTION:
% This function loads the most recently created datasets of interpolated,
% adjusted float oxygen data and interpolated glodap oxygen data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load processed float oxygen data
% obtain current day string
date_time = datetime('now');
date_time.Format = 'MMM-yyyy';
date_string = datestr(date_time,'mmm-yyyy');
% loop to subtract days until newest file is found and downloaded
count = 0; err_count = 0;
while count == err_count
    if count < 500
        try
            load(['Data/processed_float_o2_data_adjusted_' date_string float_file_ext '.mat'],...
                'float_data_adjusted','file_date');
        catch
            date_time = date_time - day(1);
            date_string = datestr(date_time,'mmm-yyyy');
            err_count = err_count + 1;
        end
        count = count + 1;
    else
        disp('Float data file could not be found.')
        return
    end
end
% display the number of matching floats and profiles
disp(['# of matching Argo profiles (OXY): ' num2str(length(unique(float_data_adjusted.OXY_PROF_ID)))]);
disp(['# of matching Argo floats (OXY): ' num2str(length(unique(float_data_adjusted.OXY_FLOAT)))]);
% clean up
clear date_time date_string count err_count

%% load processed glodap oxygen data
% obtain current day string
date_time = datetime('now');
date_time.Format = 'yyyy';
date_string = datestr(date_time,'yyyy');
% loop to subtract days until newest file is found and downloaded
count = 0; err_count = 0;
while count == err_count
    if count < 500
        try
            load(['Data/processed_glodap_o2_data_' date_string '.mat'],'glodap_data');
        catch
            date_time = date_time - day(1);
            date_string = datestr(date_time,'yyyy');
            err_count = err_count + 1;
        end
    count = count + 1;
    else
        disp('GLODAP data file could not be found.')
        return
    end
end
% display the number of matching cruises and profiles
disp(['# of matching GLODAP profiles (OXY): ' num2str(length(unique(glodap_data.OXY_ID)))]);
disp(['# of matching GLODAP cruises (OXY): ' num2str(length(unique(glodap_data.OXY_CRU)))]);
% clean up
clear date_time date_string count err_count

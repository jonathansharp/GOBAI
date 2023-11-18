% load_interpolated_combined_data_to_workspace
%
% DESCRIPTION:
% This function loads the most recently created combined dataset of
% interpolated, adjusted float oxygen data concatenated with interpolated
% glodap oxygen data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load processed combined oxygen data
% obtain current day string
date_time = datetime('now');
date_time.Format = 'MMM-yyyy';
date_string = datestr(date_time,'mmm-yyyy');
% loop to subtract days until newest file is found and downloaded
count = 0; err_count = 0;
while count == err_count
    if count < 2000
        try
            load(['Data/processed_all_o2_data_' date_string float_file_ext '.mat'],...
                'all_data','file_date');
        catch
            date_time = date_time - day(1);
            date_string = datestr(date_time,'mmm-yyyy');
            err_count = err_count + 1;
        end
        count = count + 1;
    else
        disp(['Data/processed_all_data_' date_string float_file_ext '.mat']);
        disp('Float data file could not be found.')
        return
    end
end
% display the number of matching floats and profiles
disp(['# of matching profiles: ' num2str(length(unique(all_data.id)))]);
disp(['# of matching platforms (cruises+floats): ' num2str(length(unique(all_data.platform)))]);
% clean up
clear date_time date_string count err_count

%% add spatiotemporal variables
% sine-transform longitude
all_data.lon_cos = cosd(all_data.longitude-20);
% calculate date
date = datevec(all_data.time);
% extract year
all_data.year = date(:,1);
% extract day of year
date0 = date; date0(:,2:3) = 1;
all_data.day = datenum(date) - datenum(date0);
clear date date0
% calculate sine and cosine of day
all_data.day_sin = sin((2.*pi.*all_data.day)./365.25);
all_data.day_cos = cos((2.*pi.*all_data.day)./365.25);

% %% add ancillary data
% % import chlorophyll
% 
% % calculate log10 of chlorophyll
% all_data.log10_chl = log10(all_data.chlorophyll);
% % import ETOPO2
% 
% % determine bottom depth
% all_data.bottom_depth = bottom_depth(all_data.latitude,all_data.longitude);


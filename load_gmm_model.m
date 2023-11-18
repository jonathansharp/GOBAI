% load_gmm_data_clusters
%
% DESCRIPTION:
% This function loads the most recently created Gaussian Mixture Model.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load  Gaussian Mixture Model
% obtain current day string
date_time = datetime('now');
date_time.Format = 'dd-MMM-yyyy';
date_string = datestr(date_time,'dd-mmm-yyyy');
% loop to subtract days until newest file is found and downloaded
count = 0; err_count = 0;
while count == err_count
    try
        load(['Data/all_data_clusters_' date_string '.mat'],'all_data_clusters');
    catch
        date_time = date_time - day(1);
        date_string = datestr(date_time,'dd-mmm-yyyy');
        err_count = err_count + 1;
    end
    count = count + 1;
end
% clean up
clear date_time date_string count err_count

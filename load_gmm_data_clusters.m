% load_gmm_data_clusters
%
% DESCRIPTION:
% This function loads the most recently created file of cluster assignments
% and probabilities for a combined dataset of float and glodap data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/14/2023

%% load data clusters formed from Gaussian Mixture Modelling
% obtain current day string
date_time = datetime('now');
date_time.Format = 'MMM-yyyy';
date_string = datestr(date_time,'mmm-yyyy');
% loop to subtract days until newest file is found and downloaded
count = 0; err_count = 0;
while count == err_count
    if count < 2000
        try
            load(['Data/all_data_clusters_' num2str(num_clusters) '_' ...
                date_string float_file_ext '.mat'],'all_data_clusters');
        catch
            date_time = date_time - day(1);
            date_string = datestr(date_time,'mmm-yyyy');
            err_count = err_count + 1;
        end
        count = count + 1;
    else
        disp('Cluster file could not be found.')
        return
    end
end
% clean up
clear date_time date_string count err_count

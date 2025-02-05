% reduce_ctd_data
%
% DESCRIPTION:
% Because bottle measurements and CTD sensor measurements are often taken
% on the same cast, this function identifies identical casts and removes
% the ctd data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/22/2023

%function reduce_ctd_data(glodap_year)

%% load interpolated float and glodap data
load(['O2/Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');
load(['O2/Data/processed_wod_ctd_o2_data_' num2str(glodap_year) '.mat'],...
    'wod_data');

%% define station ids and variables
stations = unique(wod_data.ID);
vars = fieldnames(wod_data);

%% check each wod ctd profile for "matching" glodap profile
for f = 1:length(stations) % for each unique station id

    %% index according to station id
    idx = wod_data.ID == stations(f);

    %% check for corresponding lat, lon, and time for glodap data
    idx_lat = abs(mean(wod_data.LAT(idx))-glodap_data.LAT) <= 0.1;
    idx_lon = abs(mean(wod_data.LON(idx))-glodap_data.LON) <= 0.1;
    idx_time = abs(mean(wod_data.TIME(idx))-glodap_data.TIME) <= 1;
    idx_all = idx_lat & idx_lon & idx_time;

    %% eliminate wod ctd data if all criteria are met
    if any(idx_all) % if all criteria are met for any station 
        for v = 1:length(vars)
            wod_data.(vars{v})(idx) = [];
        end
    end

end

%% save reduced ctd data
if ~exist([pwd '/' param1 '/Data'],'dir'); mkdir([param1 '/Data']); end
save([param1 '/Data/processed_wod_ctd_reduced_' param '_data_' year '.mat'],...
    'wod_data','-v7.3');

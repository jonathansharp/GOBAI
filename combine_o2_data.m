% combine_data
%
% DESCRIPTION:
% This function concatenates processed/adjusted float data and processed
% glodap data into one data structure.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

%% load data after implementing float data adjustment
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_float_o2_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date');
load(['Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% Combine datasets
float_vars = fieldnames(float_data_adjusted);
glodap_vars = fieldnames(glodap_data);

%% Assemble index for float oxygen data
float_idx = true(size(float_data_adjusted.OXY));
for v = 1:length(float_vars)
    float_idx(isnan(float_data_adjusted.(float_vars{v}))) = 0;
end
float_idx(float_data_adjusted.OXY_PRES<0) = 0;

%% Assemble index for glodap oxygen data
glodap_idx = true(size(glodap_data.OXY));
for v = 1:length(glodap_vars)
    glodap_idx(isnan(glodap_data.(glodap_vars{v}))) = 0;
end
glodap_idx(glodap_data.OXY_PRES<0) = 0;

%% Assemble combined dataset
all_data.platform = [float_data_adjusted.OXY_FLOAT(float_idx);...
    glodap_data.OXY_CRU(glodap_idx)];
all_data.id = [float_data_adjusted.OXY_PROF_ID(float_idx);...
    glodap_data.OXY_ID(glodap_idx)];
all_data.latitude = [float_data_adjusted.OXY_LAT(float_idx);...
    glodap_data.OXY_LAT(glodap_idx)];
all_data.longitude = [float_data_adjusted.OXY_LON(float_idx);...
    glodap_data.OXY_LON(glodap_idx)];
all_data.sigma = [float_data_adjusted.OXY_SIGMA(float_idx);...
    glodap_data.OXY_SIGMA(glodap_idx)];
all_data.pressure = [float_data_adjusted.OXY_PRES(float_idx);...
    glodap_data.OXY_PRES(glodap_idx)];
all_data.time = [float_data_adjusted.OXY_TIME(float_idx);...
    glodap_data.OXY_TIME(glodap_idx)];
all_data.temperature = [float_data_adjusted.OXY_TEMP(float_idx);...
    glodap_data.OXY_TEMP(glodap_idx)];
all_data.temperature_cns = [float_data_adjusted.OXY_CNSTEMP(float_idx);...
    glodap_data.OXY_CNSTEMP(glodap_idx)];
all_data.salinity = [float_data_adjusted.OXY_SAL(float_idx);...
    glodap_data.OXY_SAL(glodap_idx)];
all_data.salinity_abs = [float_data_adjusted.OXY_ABSSAL(float_idx);...
    glodap_data.OXY_ABSSAL(glodap_idx)];
all_data.oxygen = [float_data_adjusted.OXY(float_idx);...
    glodap_data.OXY(glodap_idx)];

% transform longitude and day of year
all_data.lon_cos = cosd(all_data.longitude-20);
all_data.day_sin = sin((2.*pi.*all_data.day)/365.25);
all_data.day_cos = cos((2.*pi.*all_data.day)/365.25);

%% save combined oxygen data
if ~exist([pwd '/Data'],'dir'); mkdir('Data'); end
save(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');
clear all_data v float_vars glodap_vars
clear float_data_adjusted float_idx glodap_data glodap_idx

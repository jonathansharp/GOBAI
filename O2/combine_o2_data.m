% combine_o2_data
%
% DESCRIPTION:
% This function concatenates processed/adjusted float data and processed
% glodap data into one data structure.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 12/1/2023

function combine_o2_data(float_file_ext,glodap_year)

%% load data after implementing float data adjustment
load(['O2/Data/processed_float_o2_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date');
load(['O2/Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
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
% all_data.sigma = [float_data_adjusted.OXY_SIGMA(float_idx);...
%     glodap_data.OXY_SIGMA(glodap_idx)];
all_data.pressure = [float_data_adjusted.OXY_PRES(float_idx);...
    glodap_data.OXY_PRES(glodap_idx)];
all_data.time = [float_data_adjusted.OXY_TIME(float_idx);...
    glodap_data.OXY_TIME(glodap_idx)];
all_data.temperature = [float_data_adjusted.OXY_TEMP(float_idx);...
    glodap_data.OXY_TEMP(glodap_idx)];
% all_data.temperature_cns = [float_data_adjusted.OXY_CNSTEMP(float_idx);...
%     glodap_data.OXY_CNSTEMP(glodap_idx)];
all_data.salinity = [float_data_adjusted.OXY_SAL(float_idx);...
    glodap_data.OXY_SAL(glodap_idx)];
% all_data.salinity_abs = [float_data_adjusted.OXY_ABSSAL(float_idx);...
%     glodap_data.OXY_ABSSAL(glodap_idx)];
all_data.oxygen = [float_data_adjusted.OXY(float_idx);...
    glodap_data.OXY(glodap_idx)];

% transform longitude and day of year
all_data.lon_cos_1 = cosd(all_data.longitude-20);
all_data.lon_cos_2 = cosd(all_data.longitude-110);
date = datevec(all_data.time);
date0 = date;
date0(:,2:3) = 0;
all_data.day = datenum(date) - datenum(date0);
all_data.day_sin = sin((2.*pi.*all_data.day)/365.25);
all_data.day_cos = cos((2.*pi.*all_data.day)/365.25);
all_data.year = date(:,1);

%% plot gridded observations
% determine bin number of each test data point on 1 degree grid
lon_edges = -180:180; lon = -179.5:179.5;
lat_edges = -90:90; lat = -89.5:89.5;
pres_edges = ([0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000])';
pres = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
[~,~,Xnum] = histcounts(all_data.longitude,lon_edges);
[~,~,Ynum] = histcounts(all_data.latitude,lat_edges);
[~,~,Znum] = histcounts(all_data.pressure,pres_edges);
% accumulate 3D grid of test data point errors
subs = [Xnum, Ynum, Znum];
idx_subs = any(subs==0,2);
sz = [length(lon),length(lat),length(pres)];
all_data.gridded_oxygen = accumarray(subs(~idx_subs,:),...
    abs(all_data.oxygen(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure('visible','off'); hold on
m_proj('robinson','lon',[20 380]);
%lon = convert_lon(lon);
[lon_temp,z] = reformat_lon(lon,all_data.gridded_oxygen(:,:,2),20);
set(gcf,'units','inches','position',[0 5 20 10]);
m_pcolor([lon_temp lon_temp(end)+1],lat,[z;z(end,:)]');
m_coast('patch',rgb('grey'));
m_grid('linestyle','-','linewidth',0.5,'xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
cmap = cmocean('ice'); cmap(1,:) = 1; colormap(cmap);
caxis([0 400]);
c=colorbar('location','southoutside');
c.Label.String = 'Average Gridded [O_{2}]';
c.FontSize = 22;
c.TickLength = 0;
if ~isfolder([pwd '/O2/Figures/Surface_Plots']); mkdir('O2/Figures/Surface_Plots'); end
exportgraphics(gcf,[pwd '/O2/Figures/Surface_Plots/Gridded_O2_10dbar.png']);
% clean up
clear land cmap c
close

%% save combined oxygen data
if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
save(['O2/Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

clear

end
% combine_data
%
% DESCRIPTION:
% This function concatenates processed/adjusted float data and processed
% glodap data into one data structure.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 12/1/2023

function combine_data(param,float_file_ext,file_date,glodap_year)

%% process parameter name
[param1,param2,param3,~,param5] = param_name(param);

%if exist([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],'file') ~= 2

%% load data after implementing float data adjustment
load([param1 '/Data/processed_float_' param '_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date');
load([param1 '/Data/processed_glodap_' param '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% Combine datasets
float_vars = fieldnames(float_data_adjusted);
glodap_vars = fieldnames(glodap_data);

%% Assemble index for float oxygen data
float_idx = true(size(float_data_adjusted.(param5)));
for v = 1:length(float_vars)
    float_idx(isnan(float_data_adjusted.(float_vars{v}))) = 0;
end
float_idx(float_data_adjusted.PRES<0) = 0;

%% Assemble index for glodap oxygen data
glodap_idx = true(size(glodap_data.(param5)));
for v = 1:length(glodap_vars)
    glodap_idx(isnan(glodap_data.(glodap_vars{v}))) = 0;
end
glodap_idx(glodap_data.PRES<0) = 0;

%% Assemble combined dataset
all_data.platform = [float_data_adjusted.(param5)(float_idx);...
    glodap_data.CRU(glodap_idx)];
all_data.id = [float_data_adjusted.PROF_ID(float_idx);...
    glodap_data.ID(glodap_idx)];
all_data.latitude = [float_data_adjusted.LAT(float_idx);...
    glodap_data.LAT(glodap_idx)];
all_data.longitude = [float_data_adjusted.LON(float_idx);...
    glodap_data.LON(glodap_idx)];
all_data.sigma = [float_data_adjusted.SIGMA(float_idx);...
    glodap_data.SIGMA(glodap_idx)];
all_data.pressure = [float_data_adjusted.PRES(float_idx);...
    glodap_data.PRES(glodap_idx)];
all_data.time = [float_data_adjusted.TIME(float_idx);...
    glodap_data.TIME(glodap_idx)];
all_data.temperature = [float_data_adjusted.TEMP(float_idx);...
    glodap_data.TEMP(glodap_idx)];
all_data.temperature_cns = [float_data_adjusted.CNSTEMP(float_idx);...
    glodap_data.CNSTEMP(glodap_idx)];
all_data.salinity = [float_data_adjusted.SAL(float_idx);...
    glodap_data.SAL(glodap_idx)];
all_data.salinity_abs = [float_data_adjusted.ABSSAL(float_idx);...
    glodap_data.ABSSAL(glodap_idx)];
all_data.(param2) = [float_data_adjusted.(param5)(float_idx);...
    glodap_data.(param5)(glodap_idx)];

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

% calculate depth
all_data.depth = -gsw_z_from_p(all_data.pressure,all_data.latitude);

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
all_data.(['gridded_' param2]) = accumarray(subs(~idx_subs,:),...
    abs(all_data.(param2)(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure('visible','off'); hold on
m_proj('robinson','lon',[20 380]);
%lon = convert_lon(lon);
[lon_temp,z] = reformat_lon(lon,all_data.(['gridded_' param2])(:,:,2),20);
set(gcf,'units','inches','position',[0 5 20 10]);
m_pcolor([lon_temp lon_temp(end)+1],lat,[z;z(end,:)]');
m_coast('patch',rgb('grey'));
m_grid('linestyle','-','linewidth',0.5,'xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
cmap = cmocean('ice'); cmap(1,:) = 1; colormap(cmap);
caxis([0 400]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Gridded ' param3];
c.FontSize = 22;
c.TickLength = 0;
if ~isfolder([pwd '/' param1 '/Figures/Surface_Plots']); mkdir([param1 '/Figures/Surface_Plots']); end
exportgraphics(gcf,[pwd '/' param1 '/Figures/Surface_Plots/Gridded_' param1 '_10dbar_' file_date float_file_ext '.png']);
% clean up
clear land cmap c
close

%% save combined oxygen data
if ~exist([pwd '/' param1 '/Data'],'dir'); mkdir([param1 '/Data']); end
save([param1 '/Data/processed_all_' param '_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

clear

%end

end
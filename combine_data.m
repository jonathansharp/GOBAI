% combine_data
%
% DESCRIPTION:
% This function concatenates processed/adjusted float data and processed
% glodap data into one data structure.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/5/2025

function combine_data(param_props,float_file_ext,glodap_year,snap_date)

%% change temporary param name for DIC
if strcmp(param_props.temp_name,'PH')
    param_props.temp_name = 'DIC';
end

%% load data after implementing float data adjustment
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load([param_props.dir_name '/Data/processed_float_' param_props.file_name '_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date');
load([param_props.dir_name '/Data/processed_glodap_' param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% Combine datasets
float_vars = fieldnames(float_data_adjusted);
glodap_vars = fieldnames(glodap_data);

%% Assemble index for float oxygen data
float_idx = true(size(float_data_adjusted.(param_props.temp_name)));
for v = 1:length(float_vars)
    float_idx(isnan(float_data_adjusted.(float_vars{v}))) = 0;
end
float_idx(float_data_adjusted.PRES<0) = 0;

%% Assemble index for glodap oxygen data
glodap_idx = true(size(glodap_data.(param_props.temp_name)));
for v = 1:length(glodap_vars)
    glodap_idx(isnan(glodap_data.(glodap_vars{v}))) = 0;
end
glodap_idx(glodap_data.PRES<0) = 0;

%% Assemble combined dataset
all_data.platform = [float_data_adjusted.(param_props.temp_name)(float_idx);...
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
all_data.(param_props.file_name) = [float_data_adjusted.(param_props.temp_name)(float_idx);...
    glodap_data.(param_props.temp_name)(glodap_idx)];
if strcmp(param_props.file_name,'dic')
    all_data.ph = [float_data_adjusted.PH(float_idx);...
        glodap_data.PH(glodap_idx)];
    all_data.talk = [float_data_adjusted.TA(float_idx);...
        glodap_data.TA(glodap_idx)];
    all_data.o2 = [float_data_adjusted.OXY(float_idx);...
        glodap_data.OXY(glodap_idx)];
    all_data.no3 = [float_data_adjusted.NIT(float_idx);...
        glodap_data.NIT(glodap_idx)];
elseif strcmp(param_props.file_name,'no3')
    all_data.o2 = [float_data_adjusted.OXY(float_idx);...
        glodap_data.OXY(glodap_idx)];
end

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
all_data.(['gridded_' param_props.file_name]) = accumarray(subs(~idx_subs,:),...
    abs(all_data.(param_props.file_name)(~idx_subs)),sz,@nanmean);
clear subs sz
% plot map
figure('visible','off'); hold on
m_proj('robinson','lon',[20 380]);
%lon = convert_lon(lon);
[lon_temp,z] = reformat_lon(lon,all_data.(['gridded_' param_props.file_name])(:,:,2),20);
set(gcf,'units','inches','position',[0 5 20 10]);
m_pcolor([lon_temp lon_temp(end)+1],lat,[z;z(end,:)]');
m_coast('patch',rgb('grey'));
m_grid('linestyle','-','linewidth',0.5,'xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
cmap = param_props.cmap; cmap(1,:) = 1; colormap(cmap);
caxis([param_props.edges(1) param_props.edges(end)]);
c=colorbar('location','southoutside');
c.Label.String = ['Average Gridded ' param_props.label];
c.FontSize = 22;
c.TickLength = 0;
if ~isfolder([pwd '/' param_props.dir_name '/Figures/Surface_Plots']); mkdir([param_props.dir_name '/Figures/Surface_Plots']); end
export_fig(gcf,[pwd '/' param_props.dir_name '/Figures/Surface_Plots/Gridded_' param_props.dir_name '_10dbar_' file_date float_file_ext '.png'],'-transparent');
% clean up
clear land cmap c
close

%% save combined oxygen data
if ~exist([pwd '/' param_props.dir_name '/Data'],'dir'); mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

end
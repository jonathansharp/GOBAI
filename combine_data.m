% combine_data
%
% DESCRIPTION:
% This function concatenates processed/adjusted float data and processed
% glodap data into one data structure.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/5/2025

function combine_data(param_props,float_file_ext,start_year,glodap_year,snap_date,...
    include_float,include_glodap,include_ctd)

%% change temporary param name for DIC
if strcmp(param_props.temp_name,'PH')
    param_props.temp_name = 'DIC';
end

%% define dataset extensions
if include_float == 1; float_ext = 'f'; else float_ext = ''; end
if include_glodap == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if include_ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% load data after implementing float data adjustment
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load([param_props.dir_name '/Data/processed_float_' param_props.file_name '_data_adjusted_' file_date float_file_ext '.mat'],...
    'float_data_adjusted','file_date');
load([param_props.dir_name '/Data/processed_glodap_' param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');
load(['O2/Data/processed_wod_ctd_' param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'wod_data');

%% Combine datasets
float_vars = fieldnames(float_data_adjusted);
glodap_vars = fieldnames(glodap_data);
wod_vars = fieldnames(wod_data);

%% Assemble index for float data
if include_float == 0
    float_idx = false(size(float_data_adjusted.(param_props.temp_name)));
else
    float_idx = true(size(float_data_adjusted.(param_props.temp_name)));
    for v = 1:length(float_vars)
        float_idx(isnan(float_data_adjusted.(float_vars{v}))) = false;
    end
    float_idx(float_data_adjusted.PRES<0) = false;
    % remove data earlier than start year
    % float_idx(float_data_adjusted.YEAR<start_year) = false;
end

%% Assemble index for glodap data
if include_glodap == 0
    glodap_idx = false(size(glodap_data.(param_props.temp_name)));
else
    glodap_idx = true(size(glodap_data.(param_props.temp_name)));
    for v = 1:length(glodap_vars)
        glodap_idx(isnan(glodap_data.(glodap_vars{v}))) = false;
    end
    glodap_idx(glodap_data.PRES<0) = false;
    % remove data earlier than start year
    % glodap_idx(glodap_data.YEAR<start_year) = false;
end

%% Assemble index for wod data
if include_ctd == 0
    wod_idx = false(size(wod_data.(param_props.temp_name)));
else
    wod_idx = true(size(wod_data.(param_props.temp_name)));
    for v = 1:length(wod_vars)
        wod_idx(isnan(wod_data.(wod_vars{v}))) = false;
    end
    wod_idx(wod_data.PRES<0) = false;
    % remove data earlier than start year
    % wod_idx(wod_data.YEAR<start_year) = false;
end

%% Assemble combined dataset
all_data.type = [repmat(1,sum(float_idx),1);...
    repmat(2,sum(glodap_idx),1);repmat(3,sum(wod_idx),1)];
all_data.platform = [float_data_adjusted.FLOAT(float_idx);...
    glodap_data.CRU(glodap_idx);wod_data.CRU(wod_idx)];
all_data.id = [float_data_adjusted.PROF_ID(float_idx);...
    glodap_data.ID(glodap_idx);wod_data.ID(wod_idx)];
all_data.latitude = [float_data_adjusted.LAT(float_idx);...
    glodap_data.LAT(glodap_idx);wod_data.LAT(wod_idx)];
all_data.longitude = [float_data_adjusted.LON(float_idx);...
    glodap_data.LON(glodap_idx);wod_data.LON(wod_idx)];
all_data.sigma = [float_data_adjusted.SIGMA(float_idx);...
    glodap_data.SIGMA(glodap_idx);wod_data.SIGMA(wod_idx)];
all_data.pressure = [float_data_adjusted.PRES(float_idx);...
    glodap_data.PRES(glodap_idx);wod_data.PRES(wod_idx)];
all_data.time = [float_data_adjusted.TIME(float_idx);...
    glodap_data.TIME(glodap_idx);wod_data.TIME(wod_idx)];
all_data.temperature = [float_data_adjusted.TEMP(float_idx);...
    glodap_data.TEMP(glodap_idx);wod_data.TEMP(wod_idx)];
all_data.temperature_cns = [float_data_adjusted.CNSTEMP(float_idx);...
    glodap_data.CNSTEMP(glodap_idx);wod_data.CNSTEMP(wod_idx)];
all_data.salinity = [float_data_adjusted.SAL(float_idx);...
    glodap_data.SAL(glodap_idx);wod_data.SAL(wod_idx)];
all_data.salinity_abs = [float_data_adjusted.ABSSAL(float_idx);...
    glodap_data.ABSSAL(glodap_idx);wod_data.ABSSAL(wod_idx)];
all_data.(param_props.file_name) = ...
    [float_data_adjusted.(param_props.temp_name)(float_idx);...
    glodap_data.(param_props.temp_name)(glodap_idx);...
    wod_data.(param_props.temp_name)(wod_idx)];
if strcmp(param_props.file_name,'dic')
    all_data.ph = [float_data_adjusted.PH(float_idx);...
        glodap_data.PH(glodap_idx);wod_data.PH(wod_idx)];
    all_data.talk = [float_data_adjusted.TA(float_idx);...
        glodap_data.TA(glodap_idx);wod_data.TA(wod_idx)];
    all_data.o2 = [float_data_adjusted.OXY(float_idx);...
        glodap_data.OXY(glodap_idx);wod_data.OXY(wod_idx)];
    all_data.no3 = [float_data_adjusted.NIT(float_idx);...
        glodap_data.NIT(glodap_idx);wod_data.NIT(wod_idx)];
elseif strcmp(param_props.file_name,'no3')
    all_data.o2 = [float_data_adjusted.OXY(float_idx);...
        glodap_data.OXY(glodap_idx);wod_data.OXY(wod_idx)];
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
% average parameter
all_data.(['gridded_' param_props.file_name]) = accumarray(subs(~idx_subs,:),...
    abs(all_data.(param_props.file_name)(~idx_subs)),sz,@nanmean);
% months represented
all_data.(['num_' param_props.file_name '_obs']) = accumarray(subs(~idx_subs,:),...
    abs(all_data.(param_props.file_name)(~idx_subs)),sz,@length);
clear subs sz
% plot map of means on each pressure level
tic; parpool(20); fprintf('Pool initiation: '); toc;
parfor lev = 1:length(pres)
    figure('visible','off'); hold on
    m_proj('robinson','lon',[20 380]);
    %lon = convert_lon(lon);
    [lon_temp,z] = reformat_lon(lon,all_data.(['gridded_' param_props.file_name])(:,:,lev),20);
    set(gcf,'units','inches','position',[0 5 20 10]);
    m_pcolor([lon_temp lon_temp(end)+1],lat,[z;z(end,:)]');
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','linewidth',0.5,'xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    cmap = param_props.cmap; cmap(1,:) = 1; colormap(cmap);
    clim([param_props.edges(1) param_props.edges(end)]);
    c=colorbar('location','southoutside');
    c.Label.String = ['Average Gridded ' param_props.label ' at ' num2str(pres(lev)) ' decibars.'];
    c.FontSize = 22;
    c.TickLength = 0;
    if ~isfolder([pwd '/' param_props.dir_name '/Figures/Surface_Plots']); mkdir([param_props.dir_name '/Figures/Surface_Plots']); end
    export_fig(gcf,[pwd '/' param_props.dir_name '/Figures/Surface_Plots/Gridded_' ...
        param_props.dir_name '_' num2str(pres(lev)) 'dbar_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.png'],'-transparent');
    close
end
% plot map of counts on each pressure level
parfor lev = 1:length(pres)
    figure('visible','off'); hold on
    m_proj('robinson','lon',[20 380]);
    % lon = convert_lon(lon);
    [lon_temp,z] = reformat_lon(lon,all_data.(['num_' param_props.file_name '_obs'])(:,:,lev),20);
    set(gcf,'units','inches','position',[0 5 20 10]);
    % set(gca,'ColorScale','log10')
    m_pcolor([lon_temp lon_temp(end)+1],lat,log10([z;z(end,:)])');
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','-','linewidth',0.5,'xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
    cmap = cmocean('amp'); cmap(1,:) = 1; colormap(cmap);
    clim([0 3]);
    c=colorbar('location','southoutside');
    c.Label.String = ['Number of ' param_props.label ' observations at ' num2str(pres(lev)) ' decibars.'];
    c.Ticks = [0 1 2 3];
    c.TickLabels = {'0' '10' '100' '1000'};
    c.FontSize = 22;
    c.TickLength = 0;
    if ~isfolder([pwd '/' param_props.dir_name '/Figures/Surface_Plots']); mkdir([param_props.dir_name '/Figures/Surface_Plots']); end
    export_fig(gcf,[pwd '/' param_props.dir_name '/Figures/Surface_Plots/' ...
        param_props.dir_name '_Obs_' num2str(pres(lev)) 'dbar_' float_ext glodap_ext ctd_ext ...
        '_' file_date float_file_ext '.png'],'-transparent');
    close
end
% end parallel session
delete(gcp('nocreate'));

%% save combined oxygen data
if ~exist([pwd '/' param_props.dir_name '/Data'],'dir'); mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' ...
    float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

end
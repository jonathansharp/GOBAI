%% Plot GOBAI versus other climatologies

%% GOBAI file information
ver = 'HR-v202606'; % version
var = 'DIC'; % variable
varname = 'dic';
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
dstn = '/fast6/dic/';
path = [dstn 'GOBAI/RFROM/FFNN/c20_Jun-2026_D/train80_val10_test10/fg/'];
fpath = [path 'GOBAI-Monthly-Mean-1x1-' var '-' ver '.nc'];
glodap_name = 'TCO2';
mobo_name = 'DIC';

%% download GLODAP gridded file
if ~isdir([dstn '/GLODAPv2.2016b_MappedClimatologies'])
    glodap_fname = 'https://www.nodc.noaa.gov/archive/arc0107/0162565/1.1/data/0-data/mapped/GLODAPv2.2016b_MappedClimatologies.tar.gz';
    websave([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz'],glodap_fname);
    gunzip([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz']);
    untar([dstn 'GLODAPv2.2016b_MappedClimatologies.tar'],dstn);
    delete([dstn 'GLODAPv2.2016b_MappedClimatologies.tar.gz']);
    delete([dstn 'GLODAPv2.2016b_MappedClimatologies.tar']);
end

%% download MOBO-DIC gridded file
if ~isfile([dstn '/MPI_MOBO-DIC_2004-2019_v2.nc'])
    mobo_fname = 'https://www.ncei.noaa.gov/data/oceans/archive/arc0211/0277099/2.3/data/0-data/MPI_MOBO-DIC_2004-2019_v2.nc';
    websave([dstn 'MPI_MOBO-DIC_2004-2019_v2.nc'],mobo_fname);
end

%% load gobai
GOBAI.lon = double(ncread(fpath,'longitude'));
GOBAI.lon = convert_lon(GOBAI.lon,'format','0-360');
GOBAI.lat = double(ncread(fpath,'latitude'));
GOBAI.pres = double(ncread(fpath,'mean_pressure'));
GOBAI.time = double(ncread(fpath,'time'));
GOBAI.var = ncread(fpath,varname);
GOBAI.var_mean = double(mean(GOBAI.var,4,'omitnan'));
GOBAI = rmfield(GOBAI,'var');

%% load glodap 
GLODAP.lon = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'lon');
GLODAP.lon = convert_lon(GLODAP.lon,'format','0-360');
GLODAP.lat = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'lat');
GLODAP.depth = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],'Depth');
GLODAP.var = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'],glodap_name);

%% load mobo-dic
mobo_local_fname = [dstn 'MPI_MOBO-DIC_2004-2019_v2.nc'];
MOBO.lon = ncread(mobo_local_fname,'lon');
MOBO.lon = convert_lon(MOBO.lon,'format','0-360');
MOBO.lat = ncread(mobo_local_fname,'lat');
MOBO.depth = ncread(mobo_local_fname,'depth');
MOBO.var = ncread(mobo_local_fname,mobo_name);
MOBO.var_mean = mean(MOBO.var,4,'omitnan');

%% calculate pressure with glodap and mobo-dic
% glodap
[GLODAP.lon3d,GLODAP.lat3d,GLODAP.depth3d] = ndgrid(GLODAP.lon,GLODAP.lat,GLODAP.depth);
% GLODAP.t = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.temperature.nc'],'temperature');
% GLODAP.s = ncread([dstn 'GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.salinity.nc'],'salinity');
GLODAP.pres3d = -gsw_p_from_z(GLODAP.depth3d,GLODAP.lat3d);
% mobo-dic
[MOBO.lon3d,MOBO.lat3d,MOBO.depth3d] = ndgrid(MOBO.lon,MOBO.lat,MOBO.depth);
% RG.t = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Temperature');
% RG.s = ncread('Data/RG_CLIM/RG_Climatology_Sal.nc','Salinity');
MOBO.pres3d = -gsw_p_from_z(MOBO.depth3d,MOBO.lat3d);

%% interpolate glodap to gobai grid
idx = ~isnan(GLODAP.var);
f = scatteredInterpolant(GLODAP.lon3d(idx),GLODAP.lat3d(idx),GLODAP.pres3d(idx),GLODAP.var(idx),'linear');
[GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d] = ndgrid(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GLODAP.var_interp = f(GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d);
GLODAP.var_interp(isnan(GOBAI.var_mean)) = NaN;

%% interpolate mean mobo-dic to gobai grid
idx = ~isnan(MOBO.var_mean);
f = scatteredInterpolant(MOBO.lon3d(idx),MOBO.lat3d(idx),MOBO.pres3d(idx),MOBO.var_mean(idx),'linear');
[GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d] = ndgrid(GOBAI.lon,GOBAI.lat,GOBAI.pres);
MOBO.var_mean_interp = f(GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d);
MOBO.var_mean_interp(isnan(GOBAI.var_mean)) = NaN;
MOBO.var_mean_interp(MOBO.var_mean_interp<0) = NaN;

%% plot GOBAI
prs = 100;
prs_idx = find(GOBAI.pres == prs);
figure('position',[100 100 800 400],'Color','w'); hold on;
title(['GOBAI-DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
worldmap([-90 90],[20 380]); mlabel off;
pcolorm(GOBAI.lat,GOBAI.lon,GOBAI.var_mean(:,:,prs_idx)');
shading flat; clim([1800 2400]); colormap(jet);
colorbar; plot_land('map');
export_fig(['Figures/GOBAI_' var '_Mean_' num2str(prs) 'dbar.png']);

%% plot GLODAP
prs = 100;
prs_idx = find(GOBAI.pres == prs);
figure('position',[100 100 800 400],'Color','w'); hold on;
title(['GLODAP DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
worldmap([-90 90],[20 380]); mlabel off;
pcolorm(GOBAI.lat,GOBAI.lon,GLODAP.var_interp(:,:,prs_idx)');
shading flat; clim([1800 2400]); colormap(jet);
colorbar; plot_land('map');
export_fig(['Figures/GLODAP_' var '_Mean_' num2str(prs) 'dbar.png']);

%% plot MOBO-DIC
prs = 100;
prs_idx = find(GOBAI.pres == prs);
figure('position',[100 100 800 400],'Color','w'); hold on;
title(['MOBO-DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
worldmap([-90 90],[20 380]); mlabel off;
pcolorm(GOBAI.lat,GOBAI.lon,MOBO.var_mean_interp(:,:,prs_idx)');
shading flat; clim([1800 2400]); colormap(jet);
colorbar; plot_land('map');
export_fig(['Figures/MOBO_' var '_Mean_' num2str(prs) 'dbar.png']);

%% plot GOBAI - GLODAP
prs = 100;
prs_idx = find(GOBAI.pres == prs);
figure('position',[100 100 800 400],'Color','w'); hold on;
title(['GOBAI - GLODAP DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
worldmap([-90 90],[20 380]); mlabel off;
pcolorm(GOBAI.lat,GOBAI.lon,GOBAI.var_mean(:,:,prs_idx)'-GLODAP.var_interp(:,:,prs_idx)');
shading flat; clim([-100 100]); colormap(cmocean('balance','pivot',0));
colorbar; plot_land('map');
export_fig(['Figures/GOBAIminusGLODAP_' var '_Mean_' num2str(prs) 'dbar.png']);

%% plot GOBAI - MOBO
prs = 100;
prs_idx = find(GOBAI.pres == prs);
figure('position',[100 100 800 400],'Color','w'); hold on;
title(['GOBAI - MOBO DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
worldmap([-90 90],[20 380]); mlabel off;
pcolorm(GOBAI.lat,GOBAI.lon,GOBAI.var_mean(:,:,prs_idx)'-MOBO.var_mean_interp(:,:,prs_idx)');
shading flat; clim([-100 100]); colormap(cmocean('balance','pivot',0));
colorbar; plot_land('map');
export_fig(['Figures/GOBAIminusMOBO_' var '_Mean_' num2str(prs) 'dbar.png']);
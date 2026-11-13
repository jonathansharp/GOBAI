%% Plot GOBAI versus other climatologies

%% GOBAI file information
ver = 'HR-v202608'; % version
var = 'DIC'; % variable
varname = 'dic';
% path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
dstn = '/fast6/dic/';
modifier = 'fg-yes_pH_corr-no_dic_adj';
path = [dstn 'GOBAI/RFROM/FFNN/c20_Aug-2026_D/train80_val10_test10/' modifier '/'];
fpath = [path 'GOBAI-Monthly-Mean-1x1-' var '-' ver '.nc'];
glodap_name = 'TCO2';
mobo_name = 'DIC';
pressures = [10 100 200 400 800 1250 1800];

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
GOBAI.year = datetime(datenum(1950,1,1)+GOBAI.time,'ConvertFrom','datenum','Format','yyyy');
GOBAI.var = ncread(fpath,varname);
time_idx_2002 = year(GOBAI.year) >= 1993 & year(GOBAI.year) <= 2011; % center on 2004
GOBAI.var_mean_2002 = double(mean(GOBAI.var(:,:,:,time_idx_2002),4,'omitnan'));
time_idx_2011 = year(GOBAI.year) >= 2004 & year(GOBAI.year) <= 2018; % center on 2011
GOBAI.var_mean_2011 = double(mean(GOBAI.var(:,:,:,time_idx_2011),4,'omitnan'));
GOBAI = rmfield(GOBAI,'var');
idx_rm_lon = find(GOBAI.lon>-5&GOBAI.lon<45);
idx_rm_lat = find(GOBAI.lat>15&GOBAI.lat<50);
GOBAI.var_mean_2002(idx_rm_lon,idx_rm_lat,:) = NaN;
GOBAI.var_mean_2011(idx_rm_lon,idx_rm_lat,:) = NaN;

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
MOBO.year = datetime(datenum(2004,ncread(mobo_local_fname,'juld'),1),'ConvertFrom','datenum','Format','yyyy');
MOBO.var = ncread(mobo_local_fname,mobo_name);
MOBO.var(MOBO.var<0) = NaN;
time_idx_2011 = year(MOBO.year) >= 2004 & year(MOBO.year) <= 2018; % center on 2011
MOBO.var_mean = mean(MOBO.var(:,:,:,time_idx_2011),4,'omitnan');

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
GLODAP.var_interp(isnan(GOBAI.var_mean_2011)) = NaN;

%% interpolate mean mobo-dic to gobai grid
idx = ~isnan(MOBO.var_mean);
f = scatteredInterpolant(MOBO.lon3d(idx),MOBO.lat3d(idx),MOBO.pres3d(idx),MOBO.var_mean(idx),'linear');
[GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d] = ndgrid(GOBAI.lon,GOBAI.lat,GOBAI.pres);
MOBO.var_mean_interp = f(GOBAI.lon3d,GOBAI.lat3d,GOBAI.pres3d);
MOBO.var_mean_interp(isnan(GOBAI.var_mean_2011)) = NaN;
MOBO.var_mean_interp(MOBO.var_mean_interp<0) = NaN;

%% plot GOBAI
for p = 1:length(pressures)
    prs = pressures(p);
    prs_idx = find(GOBAI.pres == prs);
    figure('position',[100 100 800 400],'Color','w','Visible','off'); hold on;
    set(gca,'fontsize',12);
    title(['GOBAI-DIC (2011) at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
    worldmap([-90 90],[20 380]); mlabel off; plabel off;
    pcolorm(GOBAI.lat-0.5,GOBAI.lon-0.5,GOBAI.var_mean_2011(:,:,prs_idx)');
    shading flat; clim([1900 2300]);
    colorbar; plot_land('map'); colormap(cmocean('matter'));
    export_fig(['Figures/GOBAI_' var '_Mean_' num2str(prs) 'dbar.png']);
    close;
end

%% plot GLODAP
for p = 1:length(pressures)
    prs = pressures(p);
    prs_idx = find(GOBAI.pres == prs);
    figure('position',[100 100 800 400],'Color','w','Visible','off'); hold on;
    set(gca,'fontsize',12);
    title(['GLODAP DIC (2002) at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
    worldmap([-90 90],[20 380]); mlabel off; plabel off;
    pcolorm(GOBAI.lat-0.5,GOBAI.lon-0.5,GLODAP.var_interp(:,:,prs_idx)');
    shading flat; clim([1900 2300]); colormap(jet);
    colorbar; plot_land('map'); colormap(cmocean('matter'));
    export_fig(['Figures/GLODAP_' var '_Mean_' num2str(prs) 'dbar.png']);
    close;
end

%% plot MOBO-DIC
for p = 1:length(pressures)
    prs = pressures(p);
    prs_idx = find(GOBAI.pres == prs);
    figure('position',[100 100 800 400],'Color','w','Visible','off'); hold on;
    set(gca,'fontsize',12);
    title(['MOBO-DIC (2011) at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
    worldmap([-90 90],[20 380]); mlabel off; plabel off;
    pcolorm(GOBAI.lat-0.5,GOBAI.lon-0.5,MOBO.var_mean_interp(:,:,prs_idx)');
    shading flat; clim([1900 2300]);
    colorbar; plot_land('map'); colormap(cmocean('matter'));
    export_fig(['Figures/MOBO_' var '_Mean_' num2str(prs) 'dbar.png']);
    close;
end

%% plot GOBAI - GLODAP
for p = 1:length(pressures)
    prs = pressures(p);
    prs_idx = find(GOBAI.pres == prs);
    figure('position',[100 100 800 600],'Color','w','Visible','off'); hold on;
    ax1 = worldmap([-90 90],[20 380]); mlabel off; plabel off;
    set(ax1,'fontsize',12,'Position',[0.1 0.35 0.8 0.6]);
    title(['GOBAI (2002) - GLODAP (2002) DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})']);
    diff_data = GOBAI.var_mean_2002(:,:,prs_idx)'-GLODAP.var_interp(:,:,prs_idx)';
    fprintf('Mean difference (GOBAI - GLODAP) at %d dbar: %f +/- %f\n', prs, mean(diff_data(:),'omitnan'), std(diff_data(:),'omitnan'));
    pcolorm(GOBAI.lat-0.5,GOBAI.lon-0.5,diff_data);
    shading flat; clim([-100 100]); colormap(cmocean('balance','pivot',0));
    colorbar; plot_land('map');
    ax2 = axes('Position',[0.1 0.1 0.8 0.2]); set(ax2,'fontsize',12);
    histogram(ax2,diff_data(:),'BinWidth',5);
    xlim(ax2,[-100 100]);
    xlabel(ax2,'DIC Difference (\mumol kg^{-1})');
    ylabel(ax2,'Frequency');
    export_fig(['Figures/GOBAIminusGLODAP_' var '_Mean_' num2str(prs) 'dbar.png']);
    close;
end

%% plot GOBAI - MOBO
for p = 1:length(pressures)
    prs = pressures(p);
    prs_idx = find(GOBAI.pres == prs);
    figure('position',[100 100 800 600],'Color','w','Visible','off'); hold on;
    ax1 = worldmap([-90 90],[20 380]); mlabel off; plabel off;
    set(ax1,'fontsize',12,'Position',[0.1 0.35 0.8 0.6]);
    title(['GOBAI (2011) - MOBO (2011) DIC at ' num2str(prs) ' dbar (\mumol kg^{-1})'])
    diff_data = GOBAI.var_mean_2011(:,:,prs_idx)'-MOBO.var_mean_interp(:,:,prs_idx)';
    fprintf('Mean difference (GOBAI - MOBO) at %d dbar: %f +/- %f\n', prs, mean(diff_data(:),'omitnan'), std(diff_data(:),'omitnan'));
    pcolorm(GOBAI.lat-0.5,GOBAI.lon-0.5,diff_data);
    shading flat; clim([-100 100]); colormap(cmocean('balance','pivot',0));
    colorbar; plot_land('map');
    ax2 = axes('Position',[0.1 0.1 0.8 0.2]); set(ax2,'fontsize',12);
    histogram(ax2,diff_data(:),'BinWidth',5);
    xlim(ax2,[-100 100]);
    xlabel(ax2,'DIC Difference (\mumol kg^{-1})');
    ylabel(ax2,'Frequency');
    export_fig(['Figures/GOBAIminusMOBO_' var '_Mean_' num2str(prs) 'dbar.png']);
    close;
end
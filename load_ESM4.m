function [ESM4,xdim,ydim,zdim,tdim] = load_ESM4(fpath,snap_date)

% process date
start_year = 2004;
start_month = 1;
date_str = num2str(snap_date);
end_year = str2double(date_str(1:4));
end_month = str2double(date_str(5:6));

% define paths
path1_hist = [fpath 'historical/native_grid/'];
path1_ssp = [fpath 'ssp245/native_grid/'];
path2 = '_Omon_GFDL-ESM4_';
path3 = '_r1i1p1f1_gr_';
% load time
time_hist = ncread([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc'],'time');
time_hist = daynoleap2datenum(time_hist-1,1850);
time_hist_sal = ncread([path1_hist 'so' path2 'historical' path3 '199001-200912.nc'],'time');
time_hist_sal = daynoleap2datenum(time_hist_sal-1,1850);
time_ssp = ncread([path1_ssp 'o2' path2 'ssp245' path3 '201501-210012.nc'],'time');
time_ssp = daynoleap2datenum(time_ssp-1,1850);
% create time indices
time_min = datenum([start_year start_month 15]);
time_strt = find(abs(time_hist-time_min)==min(abs(time_hist-time_min)));
time_strt_sal = find(abs(time_hist_sal-time_min)==min(abs(time_hist_sal-time_min)));
time_max = datenum([end_year end_month 15]);
time_end = find(abs(time_ssp-time_max)==min(abs(time_ssp-time_max)));
ESM4.time = single([time_hist(time_strt:end);time_ssp(1:time_end)]);
% load other 1-D model variables
ESM4.lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc'],'lon'));
ESM4.lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc'],'lat'));
ESM4.depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc'],'lev'));
% load 4-D model variables
idx_depth = find(ESM4.depth < 2100); ESM4.depth = ESM4.depth(idx_depth);
inf = ncinfo([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc']);
ESM4.oxy = single(cat(4,ncread([path1_hist 'o2' path2 'historical' path3 '185001-201412.nc'],'o2',[1 1 1 time_strt],[360 180 max(idx_depth) inf.Dimensions(4).Length-time_strt+1]),...
    ncread([path1_ssp 'o2' path2 'ssp245' path3 '201501-210012.nc'],'o2',[1 1 1 1],[360 180 max(idx_depth) time_end])));
ESM4.theta = single(cat(4,ncread([path1_hist 'thetao' path2 'historical' path3 '185001-201412.nc'],'thetao',[1 1 1 time_strt],[360 180 max(idx_depth) inf.Dimensions(4).Length-time_strt+1]),...
    ncread([path1_ssp 'thetao' path2 'ssp245' path3 '201501-210012.nc'],'thetao',[1 1 1 1],[360 180 max(idx_depth) time_end])));
inf = ncinfo([path1_hist 'so' path2 'historical' path3 '199001-200912.nc']);
ESM4.sal = single(cat(4,ncread([path1_hist 'so' path2 'historical' path3 '199001-200912.nc'],'so',[1 1 1 time_strt_sal],[360 180 max(idx_depth) inf.Dimensions(4).Length-time_strt_sal+1]),...
    ncread([path1_hist 'so' path2 'historical' path3 '201001-201412.nc'],'so',[1 1 1 1],[360 180 max(idx_depth) 60]),...
    ncread([path1_ssp 'so' path2 'ssp245' path3 '201501-203412.nc'],'so',[1 1 1 1],[360 180 max(idx_depth) time_end])));
% ESM4.chl = cat(4,ncread([path1_hist 'chl' path2 'historical' path3 '199001-200912.nc'],'chl',[1 1 1 time_strt_sal],[360 180 max(idx_depth) inf.Dimensions(4).Length-time_strt_sal+1]),...
%     ncread([path1_hist 'chl' path2 'historical' path3 '201001-201412.nc'],'chl'),...
%     ncread([path1_ssp 'chl' path2 'ssp245' path3 '201501-203412.nc'],'chl',[1 1 1 1],[360 180 max(idx_depth) time_end]));
clear path1_hist path1_ssp path2 path3 time_hist time_hist_sal time_ssp
clear inf time_min time_strt time_strt_sal time_max time_end

% limit chlorophyll to only surface values (to match observations)
% ESM4.chl = repmat(ESM4.chl(:,:,1,:),1,1,length(ESM4.depth),1);

% convert longitude to -180 to 180
ESM4.lon(ESM4.lon>180) = ESM4.lon(ESM4.lon>180) - 360;

% mask grid cells outside the RG domain
% load RG_surface_mask
% RG_surface_mask.mask = ...
%     [RG_surface_mask.mask(341:end,:);RG_surface_mask.mask(1:340,:)];
% RG_surface_mask.mask = ...
%     repmat(RG_surface_mask.mask,1,1,length(ESM4.depth),length(ESM4.time));
% ESM4.oxy(~RG_surface_mask.mask) = NaN;
% ESM4.theta(~RG_surface_mask.mask) = NaN;
% ESM4.sal(~RG_surface_mask.mask) = NaN;
% ESM4.chl(~RG_surface_mask.mask) = NaN;
% clear RG_surface_mask

% establish dimensions
xdim = length(ESM4.lon);
ydim = length(ESM4.lat);
zdim = length(ESM4.depth);
tdim = length(ESM4.time);

% calculate day
date = datevec(double(ESM4.time));
date0 = date;
date0(:,2:3) = 1;
ESM4.day = single(datenum(date) - datenum(date0));
ESM4.year = single(date(:,1));
clear date date0

% expand coordinates to 3-D variables
ESM4.lon = repmat(ESM4.lon,1,ydim,zdim);
ESM4.lat = repmat(ESM4.lat',xdim,1,zdim);
ESM4.depth = repmat(permute(ESM4.depth,[3 2 1]),xdim,ydim,1);

% calculate pressure, absolute salinity, conservative temperature,
% potential density, and distance from shore
ESM4.pres = -gsw_p_from_z(ESM4.depth,ESM4.lat);
ESM4.abs_sal = single(nan(size(ESM4.sal)));
ESM4.cns_tmp = single(nan(size(ESM4.sal)));
ESM4.tmp = single(nan(size(ESM4.sal)));
ESM4.sigma = single(nan(size(ESM4.sal)));
ESM4.dens = single(nan(size(ESM4.sal)));
for m = 1:length(ESM4.time)
    ESM4.abs_sal(:,:,:,m) = gsw_SA_from_SP(ESM4.sal(:,:,:,m),ESM4.pres,ESM4.lon,ESM4.lat);
    ESM4.cns_tmp(:,:,:,m) = gsw_CT_from_pt(ESM4.abs_sal(:,:,:,m),ESM4.theta(:,:,:,m));
    ESM4.tmp(:,:,:,m) = gsw_t_from_pt0(ESM4.abs_sal(:,:,:,m),ESM4.theta(:,:,:,m),ESM4.pres);
    ESM4.sigma(:,:,:,m) = gsw_sigma0(ESM4.abs_sal(:,:,:,m),ESM4.cns_tmp(:,:,:,m));
    ESM4.dens(:,:,:,m) = gsw_rho(ESM4.abs_sal(:,:,:,m),ESM4.cns_tmp(:,:,:,m),ESM4.pres);
end

% reduce coordinates back to 1-D variables
ESM4.lon = squeeze(ESM4.lon(:,1,1));
ESM4.lat = squeeze(ESM4.lat(1,:,1))';
ESM4.depth = squeeze(ESM4.depth(1,1,:));
ESM4.time = squeeze(ESM4.time(1,1,1,:));
ESM4.day = squeeze(ESM4.day(1,1,1,:));
ESM4.year = squeeze(ESM4.year(1,1,1,:));

% transform day by sine and cosine
ESM4.day_sin = sin((2.*pi.*ESM4.day)./365.25);
ESM4.day_cos = cos((2.*pi.*ESM4.day)./365.25);

% transform longitude by cosine:
ESM4.lon_cos1 = cosd(ESM4.lon-20);
ESM4.lon_cos2 = cosd(ESM4.lon-110);

% calculate bottom depth
ESM4.bottom_depth = ...
    single(bottom_depth(repmat(ESM4.lat',xdim,1),repmat(ESM4.lon,1,ydim)));
ESM4.bottom_depth = repmat(ESM4.bottom_depth,1,1,zdim,tdim);

% convert oxygen from mol/m^3 to umol/kg
ESM4.oxy = (ESM4.oxy./ESM4.dens).*(10^6);

% zero trap for oxygen
ESM4.oxy(ESM4.oxy<0) = 0;

% calculate oxygen solubility
ESM4.oxy_sat = o2satv2b(ESM4.sal,ESM4.tmp);

% clean up
ESM4 = rmfield(ESM4,{'day' 'theta'});

% save processed ESM4 data
save(['Data/ESM4_processed_' date_str],'ESM4','-v7.3')

% Plot T at 20 m
figure; worldmap([-90 90],[20 380]);
title('Annual mean at 20 m (GFDL-ESM4)','fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(ESM4.lat),double([ESM4.lon;ESM4.lon(end)+1]),...
    mean([ESM4.tmp(:,:,3,:);ESM4.tmp(end,:,3,:)],4,'omitnan')');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([0 30]);
colormap(cmocean('thermal',12));
c.Label.String = ['Temperature ' char(176) 'C'];
c.FontSize = 12;
c.TickLength = 0;
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Surface_Plots','dir'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,'Figures/Surface_Plots/temp_20_m_ESM4.png');
close

% Plot S at 20 m
figure; worldmap([-90 90],[20 380]);
title('Annual mean at 20 m (GFDL-ESM4)','fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(ESM4.lat),double([ESM4.lon;ESM4.lon(end)+1]),...
    mean([ESM4.sal(:,:,3,:);ESM4.tmp(end,:,3,:)],4,'omitnan')');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([32 38]);
colormap(cmocean('haline',12));
c.Label.String = 'Salinity';
c.FontSize = 12;
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Surface_Plots','dir'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,'Figures/Surface_Plots/sal_20_m_ESM4.png');
close

% Plot O2 at 20 m
figure; worldmap([-90 90],[20 380]);
title('Annual mean at 20 m (GFDL-ESM4)','fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(ESM4.lat),double([ESM4.lon;ESM4.lon(end)+1]),...
    mean([ESM4.oxy(:,:,3,:);ESM4.oxy(end,:,3,:)],4,'omitnan')');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([150 350]);
colormap(cmocean('ice',14));
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
c.FontSize = 12;
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Surface_Plots','dir'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,'Figures/Surface_Plots/oxy_20_m_ESM4.png');
close

end

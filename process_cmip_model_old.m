function process_cmip_model_old(model,fpath,snap_date,start_year,grid_type)

% process date
start_month = 1;
date_str = num2str(snap_date);
end_year = str2double(date_str(1:4));
end_month = str2double(date_str(5:6));

% define paths
path1_hist = [fpath 'historical/native_grid/'];
path1_ssp = [fpath 'ssp245/native_grid/'];
path2 = ['_Omon_' model '_'];
path3 = ['_r1i1p1f1_' grid_type];
if ~isfolder([fpath 'combined/native_grid/']); mkdir([fpath 'combined/native_grid/']); end

% load time
time_inf_historical = ncinfo([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'time');
origin_idx = find(strcmp({time_inf_historical.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_historical.Attributes(origin_idx).Value,'days since ')));
time_hist = ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'time');
time_hist = daynoleap2datenum(time_hist-1,origin_date(1));
time_inf_ssp = ncinfo([path1_ssp 'o2' path2 'ssp245' path3 '_201501-213412.nc'],'time');
origin_idx = find(strcmp({time_inf_ssp.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_ssp.Attributes(origin_idx).Value,'days since ')));
time_ssp = ncread([path1_ssp 'o2' path2 'ssp245' path3 '_201501-213412.nc'],'time');
time_ssp = daynoleap2datenum(time_ssp-1,origin_date(1));

% create time indices
time_min = datenum([start_year start_month 15]);
time_strt = find(abs(time_hist-time_min)==min(abs(time_hist-time_min)));
% time_strt_sal = find(abs(time_hist_sal-time_min)==min(abs(time_hist_sal-time_min)));
time_max = datenum([end_year end_month 15]);
time_end = find(abs(time_ssp-time_max)==min(abs(time_ssp-time_max)));
time = single([time_hist(time_strt:end);time_ssp(1:time_end)]);

% load other 1-D model variables
if strcmp(grid_type,'gr')
    lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'lon'));
    lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'lat'));
    depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'lev'));
elseif strcmp(grid_type,'gn')
    lon = (0.5:359.5)';
    lat = (-89.5:89.5)';
    nav_lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'nav_lon'));
    nav_lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'nav_lat'));
    depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_199001-200912.nc'],'olevel'));
end

% index to above 2100m
idx_depth = find(depth < 2100); depth = depth(idx_depth);

% load, combine, and save potential temperature
if isfile([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    dims_hist = ncinfo([path1_hist 'thetao' path2 'historical' path3 '_199001-200912.nc'],'thetao');
    for t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist.Size(4)-time_strt
            thetao = ncread([path1_hist 'thetao' path2 'historical' path3 '_199001-200912.nc'],...
                'thetao',[1 1 1 time_strt+t],[dims_hist.Size(1) dims_hist.Size(2) max(idx_depth) 1]);
        else
            thetao = ncread([path1_ssp 'thetao' path2 'ssp245' path3 '_201501-210012.nc'],...
                'thetao',[1 1 1 1],[dims_hist.Size(1) dims_hist.Size(2) max(idx_depth) time_end]);
        end
        % save variable to file
        if strcmp(grid_type,'gr')
            if t ==1
                ncsave_4d([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],...
                    {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
                    {'thetao' thetao 'Sea Water Potential Temperature' 'degC'});
            else
                ncwrite([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],...
                    'thetao',thetao);
            end
        elseif strcmp(grid_type,'gn')
            % interpolate
            interpolant = scatteredInterpolant(repmat(nav_lon,1,1,length(depth)),repmat(nav_lat,1,1,length(depth)),...
                repmat(permute(depth,[3 2 1]),size(nav_lon,1),size(nav_lon,2),1),thetao);
            if t ==1    
                ncsave_4d([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],...
                    {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
                    {'thetao' thetao 'Sea Water Potential Temperature' 'degC'});
            else
                ncwrite([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],...
                    'thetao',thetao);
            end
        end
    end
end

% load, combine, and save salinity
if isfile([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    inf = ncinfo([path1_hist 'so' path2 'historical' path3 '_199001-200912.nc']);
    so = single(cat(4,ncread([path1_hist 'so' path2 'historical' path3 '_199001-200912.nc'],'so',[1 1 1 time_strt_sal],[360 180 max(idx_depth) Inf]),...
        ncread([path1_hist 'so' path2 'historical' path3 '_201001-201412.nc'],'so',[1 1 1 1],[360 180 max(idx_depth) 60]),...
        ncread([path1_ssp 'so' path2 'ssp245' path3 '_201501-203412.nc'],'so',[1 1 1 1],[360 180 max(idx_depth) time_end])));
    ncsave_4d([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'so' so 'Sea Water Salinity' ''});
    clear so
end

% % load, combine, and save chlorophyll
% chl = cat(4,ncread([path1_hist 'chl' path2 'historical' path3 '_199001-200912.nc'],'chl',[1 1 1 time_strt_sal],[360 180 max(idx_depth) Inf]),...
%     ncread([path1_hist 'chl' path2 'historical' path3 '_201001-201412.nc'],'chl'),...
%     ncread([path1_ssp 'chl' path2 'ssp245' path3 '_201501-203412.nc'],'chl',[1 1 1 1],[360 180 max(idx_depth) time_end]));

% limit chlorophyll to only surface values (to match observations)
% chl = repmat(chl(:,:,1,:),1,1,length(depth),1);

% convert longitude to -180 to 180
% lon(lon>180) = lon(lon>180) - 360;

% mask grid cells outside the RG domain
% load RG_surface_mask
% RG_surface_mask.mask = ...
%     [RG_surface_mask.mask(341:end,:);RG_surface_mask.mask(1:340,:)];
% RG_surface_mask.mask = ...
%     repmat(RG_surface_mask.mask,1,1,length(depth),length(time));
% oxy(~RG_surface_mask.mask) = NaN;
% theta(~RG_surface_mask.mask) = NaN;
% sal(~RG_surface_mask.mask) = NaN;
% chl(~RG_surface_mask.mask) = NaN;
% clear RG_surface_mask

% establish dimensions
xdim = length(lon); ydim = length(lat); zdim = length(depth);

% expand coordinates to 3-D variables
lon_3d = repmat(lon,1,ydim,zdim);
lat_3d = repmat(lat',xdim,1,zdim);
depth_3d = repmat(permute(depth,[3 2 1]),xdim,ydim,1);

% calculate 3d pressure
pres_3d = -gsw_p_from_z(depth_3d,lat_3d);

% calculate and save absolute salinity
if isfile([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    so = ncread([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'],'so');
    abs_sal = single(nan(size(so)));
    for m = 1:length(time)
        abs_sal(:,:,:,m) = gsw_SA_from_SP(so(:,:,:,m),pres_3d,lon_3d,lat_3d);
    end
    ncsave_4d([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'abs_sal' abs_sal 'Absolute Salinity' ''});
    clear so abs_sal
end

% calculate and save conservative temperature
if isfile([fpath 'combined/native_grid/cns_tmp' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/cns_tmp' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],'abs_sal');
    cns_tmp = single(nan(size(thetao)));
    for m = 1:length(time)
        cns_tmp(:,:,:,m) = gsw_CT_from_pt(abs_sal(:,:,:,m),thetao(:,:,:,m));
    end
    ncsave_4d([fpath 'combined/native_grid/cns_tmp' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'cns_tmp' cns_tmp 'Conservative Temperature' 'degC'});
    clear thetao cns_tmp
end

% calculate and save in situ temperature
if isfile([fpath 'combined/native_grid/tmp' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/tmp' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/native_grid/thetao' path2 'combined' path3 '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],'abs_sal');
    tmp = single(nan(size(thetao)));
    for m = 1:length(time)
        tmp(:,:,:,m) = gsw_t_from_pt0(abs_sal(:,:,:,m),thetao(:,:,:,m),pres_3d);
    end
    ncsave_4d([fpath 'combined/native_grid/tmp' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'tmp' tmp 'In Situ Temperature' 'degC'});
    clear thetao abs_sal tmp
end

% calculate and save potential density anomaly
if isfile([fpath 'combined/native_grid/sigma' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/sigma' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/native_grid/cns_tmp' path2 'combined' path3 '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],'abs_sal');
    sigma = single(nan(size(cns_tmp)));
    for m = 1:length(time)
        sigma(:,:,:,m) = gsw_sigma0(abs_sal(:,:,:,m),cns_tmp(:,:,:,m));
    end
    ncsave_4d([fpath 'combined/native_grid/sigma' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'sigma' sigma 'Potential Density Anomaly' 'kg/m^3'});
    clear cns_tmp abs_sal sigma
end

% calculate and save in situ density
if isfile([fpath 'combined/native_grid/dens' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/dens' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/native_grid/cns_tmp' path2 'combined' path3 '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/native_grid/abs_sal' path2 'combined' path3 '.nc'],'abs_sal');
    dens = single(nan(size(cns_tmp)));
    for m = 1:length(time)
        dens(:,:,:,m) = gsw_rho(abs_sal(:,:,:,m),cns_tmp(:,:,:,m),pres_3d);
    end
    ncsave_4d([fpath 'combined/native_grid/dens' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'dens' dens 'In Situ Density' 'kg/m^3'});
    clear cns_tmp abs_sal dens
end

% calculate and save oxygen saturation
if isfile([fpath 'combined/native_grid/o2_sat' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/o2_sat' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    tmp = ncread([fpath 'combined/native_grid/tmp' path2 'combined' path3 '.nc'],'tmp');
    so = ncread([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'],'so');
    o2_sat = single(nan(size(tmp)));
    for m = 1:length(time)
        o2_sat(:,:,:,m) = o2satv2b(so(:,:,:,m),tmp(:,:,:,m));
    end
    ncsave_4d([fpath 'combined/native_grid/o2_sat' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'o2_sat' o2_sat 'Oxygen Saturation' 'micromoles per kilogram'});
    clear tmp so o2_sat
end

% load, combine, and save o2
if isfile([fpath 'combined/native_grid/o2' path2 'combined' path3 '.nc'])
    l = nc_dim_length([fpath 'combined/native_grid/o2' path2 'combined' path3 '.nc'],'time');
else; l = 0; end
if l ~= length(time)
    inf = ncinfo([path1_hist 'o2' path2 'historical' path3 '_185001-201412.nc']);
    o2 = single(cat(4,ncread([path1_hist 'o2' path2 'historical' path3 '_185001-201412.nc'],'o2',[1 1 1 time_strt],[360 180 max(idx_depth) Inf]),...
        ncread([path1_ssp 'o2' path2 'ssp245' path3 '_201501-210012.nc'],'o2',[1 1 1 1],[360 180 max(idx_depth) time_end])));
    dens = ncread([fpath 'combined/native_grid/dens' path2 'combined' path3 '.nc'],'dens');
    o2 = (o2./dens).*(10^6); % convert oxygen from mol/m^3 to umol/kg
    ncsave_4d([fpath 'combined/native_grid/o2' path2 'combined' path3 '.nc'],...
        {'lon' lon 'longitude' 'degrees east'},{'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},{'time' time 'time' 'days since 0000-01-01'},...
        {'o2' o2 'dissolved oxygen content' 'micromol per kg'});
    clear dens o2
end

% add 3d pressure to NetCDFs
vars = {'o2' 'thetao' 'so' 'abs_sal' 'cns_tmp' 'tmp' 'sigma' 'dens'};
for v = 1:length(vars)
    if isfile([fpath 'combined/native_grid/' vars{v} path2 'combined' path3 '.nc'])
        if ~nc_var_exist([fpath 'combined/native_grid/' vars{v} path2 'combined' path3 '.nc'],'pres')
            nccreate([fpath 'combined/native_grid/' vars{v} path2 'combined' path3 '.nc'],...
                'pres','Dimensions',{'lon' length(lon) 'lat' length(lat) 'depth' length(depth)});
            ncwrite([fpath 'combined/native_grid/' vars{v} path2 'combined' path3 '.nc'],...
                'pres',pres_3d);
        end
    end
end

% Plot T at 20 m
idx_depth = find(depth == 20);
tmp = ncread([fpath 'combined/native_grid/tmp' path2 'combined' path3 '.nc'],...
    'tmp',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_tmp = mean(tmp,4,'omitnan');
figure; worldmap([-90 90],[20 380]);
title(['Annual mean at 20 m (' model ')'],'fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(lat),double([lon;lon(end)+1]),...
    [mean_tmp;mean_tmp(end,:)]');
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
exportgraphics(gcf,['Figures/Surface_Plots/temp_20m_' model '.png']);
close

% Plot S at 20 m
idx_depth = find(depth == 20);
so = ncread([fpath 'combined/native_grid/so' path2 'combined' path3 '.nc'],...
    'so',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_so = mean(so,4,'omitnan');
figure; worldmap([-90 90],[20 380]);
title(['Annual mean at 20 m (' model ')'],'fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(lat),double([lon;lon(end)+1]),...
    [mean_so;mean_so(end,:)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([32 38]);
colormap(cmocean('haline',12));
c.Label.String = 'Salinity';
c.FontSize = 12;
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Surface_Plots','dir'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,['Figures/Surface_Plots/sal_20m_' model '.png']);
close

% Plot O2 at 20 m
idx_depth = find(depth == 20);
o2 = ncread([fpath 'combined/native_grid/o2' path2 'combined' path3 '.nc'],...
    'o2',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_o2 = mean(o2,4,'omitnan');
figure; worldmap([-90 90],[20 380]);
title(['Annual mean at 20 m (' model ')'],'fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(lat),double([lon;lon(end)+1]),...
    [mean_o2;mean_o2(end,:)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([150 350]);
colormap(cmocean('ice',14));
c.Label.String = '[O_{2}] (\mumol kg^{-1})';
c.FontSize = 12;
mlabel off; plabel off;
if ~exist('Figures','dir'); mkdir('Figures'); end
if ~exist('Figures/Surface_Plots','dir'); mkdir('Figures/Surface_Plots'); end
exportgraphics(gcf,['Figures/Surface_Plots/oxy_20m_' model '.png']);
close

end

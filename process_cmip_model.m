function process_cmip_model(model,fpath,snap_date,start_year,grid_label,grid_type)

%% process date
start_month = 1;
date_str = num2str(snap_date);
end_year = str2double(date_str(1:4));
end_month = str2double(date_str(5:6));

%% define paths
path1_hist = [fpath 'historical/' grid_type '/'];
path1_ssp = [fpath 'ssp245/' grid_type '/'];
path2 = ['_Omon_' model '_'];
path3 = ['_r1i1p1f1_' grid_label];
if ~isfolder([fpath 'combined/' grid_type '/']); mkdir([fpath 'combined/' grid_type '/']); end
% time extensions for cmip model paths
if strcmp(model,'GFDL-ESM4')
    ext_hist1 = '199001-200912';
    ext_ssp1 = '201501-203412';
elseif strcmp(model,'NorESM2-LM')
    ext_hist1 = '200001-200912';
    ext_ssp1 = '201501-202012';
    ext_ssp2 = '202101-203012';
end
ext_hist2 = '201001-201412';

%% load and process historical time
time_inf_historical = ncinfo([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'time');
origin_idx = find(strcmp({time_inf_historical.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_historical.Attributes(origin_idx).Value,'days since ')));
time_hist1 = ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'time');
time_hist2 = ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist2 '.nc'],'time');
time_hist = daynoleap2datenum([time_hist1;time_hist2]-1,origin_date(1));

%% load and process ssp time
time_inf_ssp = ncinfo([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp1 '.nc'],'time');
origin_idx = find(strcmp({time_inf_ssp.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_ssp.Attributes(origin_idx).Value,'days since ')));
time_ssp1 = ncread([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp1 '.nc'],'time');
if strcmp(model,'NorESM2-LM')
    time_ssp2 = ncread([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp2 '.nc'],'time');
else; time_ssp2 = []; end
time_ssp = daynoleap2datenum([time_ssp1;time_ssp2]-1,origin_date(1));

%% create time indices
time_min = datenum([start_year start_month 15]);
time_strt = find(abs(time_hist-time_min)==min(abs(time_hist-time_min)));
time_max = datenum([end_year end_month 15]);
time_end = find(abs(time_ssp-time_max)==min(abs(time_ssp-time_max)));
time = single([time_hist(time_strt:end);time_ssp(1:time_end)]);

%% load other 1-D model variables
if strcmp(grid_label,'gr')
    lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lon'));
    lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lat'));
    depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lev'));
elseif strcmp(grid_label,'gn')
    lon = (0.5:359.5)';
    lat = (-89.5:89.5)';
    nav_lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'nav_lon'));
    nav_lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'nav_lat'));
    depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'olevel'));
end

% index to above 2100m
idx_depth = find(depth < 2100); depth = depth(idx_depth);

%% load, combine, and save potential temperature
nc_filepath = [fpath 'combined/' grid_type '/thetao' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,path1_hist,path1_ssp,path2,path3,lon,lat,...
        depth,time,time_strt,grid_label,'thetao',idx_depth,'Sea Water Potential Temperature','degC');
    % plot global mean timeseries
    thetao = ncread(nc_filepath,'thetao');
    plot_global_timeseries(lat,lon,depth,time,thetao,'thetao',...
        'Potential Temperature',[char(176) 'C'],model,fpath,grid_type);
    clear thetao
end

%% load, combine, and save salinity
nc_filepath = [fpath 'combined/' grid_type '/so' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,path1_hist,path1_ssp,path2,path3,lon,lat,...
        depth,time,time_strt,grid_label,'so',idx_depth,'Sea Water Salinity','n/a');
    % plot global mean timeseries
    so = ncread(nc_filepath,'so');
    plot_global_timeseries(lat,lon,depth,time,so,'so','Salinity',...
        'n/a',model,fpath,grid_type);
    clear so
end

%% load, combine, and save chlorophyll
% nc_filepath = [fpath 'combined/' grid_type '/chl' path2 ... % define filepath
%     'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% if l ~= length(time)
%     combine_cmip_field(nc_filepath,path1_hist,path1_ssp,lon,lat,...
%         depth,time,time_strt,grid_label,'chl',idx_depth);
%     limit chlorophyll to only surface values (to match observations)
%     chl = repmat(chl(:,:,1,:),1,1,length(depth),1);
%     % plot global mean timeseries
%     chl = ncread(nc_filepath,'chl');
%     plot_global_timeseries(lat,lon,depth,time,chl,'chl','Chlorophyll','mg/m2',model,fpath,grid_type);
% end

%% process dimensional variables
% establish dimensions
xdim = length(lon); ydim = length(lat); zdim = length(depth);
% expand coordinates to 3-D variables
lon_3d = repmat(lon,1,ydim,zdim);
lat_3d = repmat(lat',xdim,1,zdim);
depth_3d = repmat(permute(depth,[3 2 1]),xdim,ydim,1);
% calculate 3d pressure
pres_3d = -gsw_p_from_z(depth_3d,lat_3d);

%% calculate and save absolute salinity
nc_filepath = [fpath 'combined/' grid_type '/abs_sal' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    so = ncread([fpath 'combined/' grid_type '/so' path2 'combined' ...
        path3 '_' num2str(start_year) '01-' date_str '.nc'],'so');
    abs_sal = single(nan(size(so)));
    for m = 1:length(time)
        abs_sal(:,:,:,m) = gsw_SA_from_SP(so(:,:,:,m),pres_3d,lon_3d,lat_3d);
    end
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'abs_sal' abs_sal 'Absolute Salinity' ''});
    clear so abs_sal
    % plot global mean timeseries
    abs_sal = ncread(nc_filepath,'abs_sal');
    plot_global_timeseries(lat,lon,depth,time,abs_sal,'abs_sal',...
        'Absolute Salinity','n/a',model,fpath,grid_type);
    clear abs_sal
end

%% calculate and save conservative temperature
nc_filepath = [fpath 'combined/' grid_type '/cns_tmp' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/' grid_type '/thetao' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
    cns_tmp = single(nan(size(thetao)));
    for m = 1:length(time)
        cns_tmp(:,:,:,m) = gsw_CT_from_pt(abs_sal(:,:,:,m),thetao(:,:,:,m));
    end
    ncsave_4d(nc_filepath,...
        {'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'cns_tmp' cns_tmp 'Conservative Temperature' 'degC'});
    clear thetao cns_tmp
    % plot global mean timeseries
    cns_tmp = ncread(nc_filepath,'cns_tmp');
    plot_global_timeseries(lat,lon,depth,time,cns_tmp,'cns_tmp',...
        'Conservative Temperature',[char(176) 'C'],model,fpath,grid_type);
    clear cns_tmp
end

%% calculate and save in situ temperature
nc_filepath = [fpath 'combined/' grid_type '/tmp' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/' grid_type '/thetao' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
    tmp = single(nan(size(thetao)));
    for m = 1:length(time)
        tmp(:,:,:,m) = gsw_t_from_pt0(abs_sal(:,:,:,m),thetao(:,:,:,m),pres_3d);
    end
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'tmp' tmp 'In Situ Temperature' 'degC'});
    clear thetao abs_sal tmp
    % plot global mean timeseries
    tmp = ncread(nc_filepath,'tmp');
    plot_global_timeseries(lat,lon,depth,time,tmp,'tmp',...
        'In Situ Temperature',[char(176) 'C'],model,fpath,grid_type);
    clear tmp
end

%% calculate and save potential density anomaly
nc_filepath = [fpath 'combined/' grid_type '/sigma' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/' grid_type '/cns_tmp' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
    sigma = single(nan(size(cns_tmp)));
    for m = 1:length(time)
        sigma(:,:,:,m) = gsw_sigma0(abs_sal(:,:,:,m),cns_tmp(:,:,:,m));
    end
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'sigma' sigma 'Potential Density Anomaly' 'kg/m^3'});
    clear cns_tmp abs_sal sigma
    % plot global mean timeseries
    sigma = ncread(nc_filepath,'sigma');
    plot_global_timeseries(lat,lon,depth,time,sigma,'sigma',...
        'Potential Density','kg/m^3',model,fpath,grid_type);
    clear sigma
end

%% calculate and save in situ density
nc_filepath = [fpath 'combined/' grid_type '/dens' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/' grid_type '/cns_tmp' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
    dens = single(nan(size(cns_tmp)));
    for m = 1:length(time)
        dens(:,:,:,m) = gsw_rho(abs_sal(:,:,:,m),cns_tmp(:,:,:,m),pres_3d);
    end
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'dens' dens 'In Situ Density' 'kg/m^3'});
    clear cns_tmp abs_sal dens
    % plot global mean timeseries
    dens = ncread(nc_filepath,'dens');
    plot_global_timeseries(lat,lon,depth,time,dens,'dens',...
        'Density','kg/m^3',model,fpath,grid_type);
    clear dens
end

%% calculate and save oxygen saturation
nc_filepath = [fpath 'combined/' grid_type '/o2_sat' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    tmp = ncread([fpath 'combined/' grid_type '/tmp' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'tmp');
    so = ncread([fpath 'combined/' grid_type '/so' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'so');
    o2_sat = single(nan(size(tmp)));
    for m = 1:length(time)
        o2_sat(:,:,:,m) = o2satv2b(so(:,:,:,m),tmp(:,:,:,m));
    end
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'o2_sat' o2_sat 'Oxygen Saturation' 'micromoles per kilogram'});
    clear tmp so o2_sat
    % plot global mean timeseries
    o2_sat = ncread(nc_filepath,'o2_sat');
    plot_global_timeseries(lat,lon,depth,time,o2_sat,'o2_sat',...
        'Oxygen Saturation','umol/kg',model,fpath,grid_type);
    clear o2_sat
end

%% load, combine, and save o2
nc_filepath = [fpath 'combined/' grid_type '/o2' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,path1_hist,path1_ssp,path2,path3,lon,lat,...
        depth,time,time_strt,grid_label,'o2',idx_depth,'Dissolved Oxygen Content','mL/L');
    o2 = ncread(nc_filepath,'o2');
    dens = ncread([fpath 'combined/' grid_type '/dens' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'],'dens');
    o2 = (o2./dens).*(10^6); % convert oxygen from mol/m^3 to umol/kg
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'o2' o2 'Dissolved Oxygen Content','umol/kg'});
    % plot global mean timeseries
    o2 = ncread(nc_filepath,'o2');
    plot_global_timeseries(lat,lon,depth,time,o2,'o2',...
        'Oxygen Amount Content','umol/kg',model,fpath,grid_type);
    clear o2
end

%% add 3d pressure to NetCDFs
vars = {'o2' 'thetao' 'so' 'abs_sal' 'cns_tmp' 'tmp' 'sigma' 'dens'};
for v = 1:length(vars)
    if isfile([fpath 'combined/' grid_type '/' vars{v} path2 'combined' path3 '.nc'])
        if ~nc_var_exist([fpath 'combined/' grid_type '/' vars{v} path2 'combined' path3 '.nc'],'pres')
            nccreate([fpath 'combined/' grid_type '/' vars{v} path2 'combined' path3 '.nc'],...
                'pres','Dimensions',{'lon' length(lon) 'lat' length(lat) 'depth' length(depth)});
            ncwrite([fpath 'combined/' grid_type '/' vars{v} path2 'combined' path3 '.nc'],...
                'pres',pres_3d);
        end
    end
end

%% Plot T at 20 m
nc_filepath = [fpath 'combined/' grid_type '/tmp' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(depth == 20);
tmp = ncread(nc_filepath,'tmp',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_tmp = mean(tmp,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
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
if ~isfolder(['Figures/' model]); mkdir(['Figures/' model]); end
export_fig(['Figures/' model '/temp_20m.png'],'-transparent');
close

%% Plot S at 20 m
nc_filepath = [fpath 'combined/' grid_type '/so' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(depth == 20);
so = ncread(nc_filepath,'so',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_so = mean(so,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
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
if ~isfolder(['Figures/' model]); mkdir(['Figures/' model]); end
export_fig(['Figures/' model '/sal_20m.png'],'-transparent');
close

%% Plot O2 at 20 m
nc_filepath = [fpath 'combined/' grid_type '/o2' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(depth == 20);
o2 = ncread(nc_filepath,'o2',[1 1 idx_depth 1],[Inf Inf 1 Inf]);
mean_o2 = mean(o2,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
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
if ~isfolder(['Figures/' model]); mkdir(['Figures/' model]); end
export_fig(['Figures/' model '/oxy_20m.png'],'-transparent');
close

end


%% embedded functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load and combine cmip model output
function combine_cmip_field(model,nc_filepath,path1_hist,path1_ssp,path2,path3,...
    lon,lat,depth,time,time_strt,grid_label,var_name,idx_depth,long_name,units)

if strcmp(model,'GFDL-ESM4')

    dims_hist1 = ncinfo([path1_hist var_name path2 'historical' path3 '_199001-200912.nc'],var_name);
    dims_hist2 = ncinfo([path1_hist var_name path2 'historical' path3 '_201001-201412.nc'],var_name);
    
    for t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist1.Size(4)-time_strt
            var = ncread([path1_hist var_name path2 'historical' path3 '_199001-200912.nc'],...
                var_name,[1 1 1 (time_strt-1)+t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)
            var = ncread([path1_hist var_name path2 'historical' path3 '_201001-201412.nc'],...
                var_name,[1 1 1 time_strt+(t-dims_hist1.Size(4))],[dims_hist2.Size(1) dims_hist2.Size(2) max(idx_depth) 1]);
        else
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_201501-203412.nc'],...
                var_name,[1 1 1 time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4))],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        end
        % save variable to file
        if strcmp(grid_label,'gr')
            if t == 1
                % create combined file
                ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
                    {'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},...
                    {'time' time(t) 'time' 'time'},...
                    {var_name var long_name units});
            else
                % append to combined file
                ncwrite(nc_filepath,'time',time(2),t);
                ncwrite(nc_filepath,var_name,var,[1 1 1 t]);
            end
        elseif strcmp(grid_label,'gn')
            % interpolate
    %         interpolant = scatteredInterpolant(repmat(nav_lon,1,1,length(depth)),repmat(nav_lat,1,1,length(depth)),...
    %             repmat(permute(depth,[3 2 1]),size(nav_lon,1),size(nav_lon,2),1),var);
            if t == 1    
                ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
                    {'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},...
                    {'time' time 'time' 'days since 0000-01-01'},...
                    {var_name var 'Sea Water Potential Temperature' 'degC'});
            else
                ncwrite(nc_filepath,'time',time(2),t);
                ncwrite(nc_filepath,var_name,var,[1 1 1 t]);
            end
        end
    end

elseif strcmp(model,'NorESM2-LM')

    dims_hist1 = ncinfo([path1_hist var_name path2 'historical' path3 '200001-200912.nc'],var_name);
    dims_hist2 = ncinfo([path1_hist var_name path2 'historical' path3 '_201001-201412.nc'],var_name);
    dims_ssp1 = ncinfo([path1_ssp var_name path2 'ssp245' path3 '_201501-202012.nc'],var_name);

    for t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist1.Size(4)-time_strt
            var = ncread([path1_hist var_name path2 'historical' path3 '_199001-200912.nc'],...
                var_name,[1 1 1 (time_strt-1)+t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)
            var = ncread([path1_hist var_name path2 'historical' path3 '_201001-201412.nc'],...
                var_name,[1 1 1 time_strt+(t-dims_hist1.Size(4))],[dims_hist2.Size(1) dims_hist2.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)+dims_ssp1.Size(4)
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_201501-202012.nc'],...
                var_name,[1 1 1 time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4))],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        else
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_202101-203012.nc'],...
                var_name,[1 1 1 time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4)-dims_ssp1.Size(4))],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        end
        % save variable to file
        if strcmp(grid_label,'gr')
            if t == 1
                % create combined file
                ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
                    {'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},...
                    {'time' time(t) 'time' 'time'},...
                    {var_name var long_name units});
            else
                % append to combined file
                ncwrite(nc_filepath,'time',time(2),t);
                ncwrite(nc_filepath,var_name,var,[1 1 1 t]);
            end
        elseif strcmp(grid_label,'gn')
            % interpolate
    %         interpolant = scatteredInterpolant(repmat(nav_lon,1,1,length(depth)),repmat(nav_lat,1,1,length(depth)),...
    %             repmat(permute(depth,[3 2 1]),size(nav_lon,1),size(nav_lon,2),1),var);
            if t == 1    
                ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
                    {'lat' lat 'latitude' 'degrees north'},...
                    {'depth' depth 'depth' 'meters'},...
                    {'time' time 'time' 'days since 0000-01-01'},...
                    {var_name var 'Sea Water Potential Temperature' 'degC'});
            else
                ncwrite(nc_filepath,'time',time(2),t);
                ncwrite(nc_filepath,var_name,var,[1 1 1 t]);
            end
        end
    end

end

end

%% plot global timeseries
function plot_global_timeseries(lat,lon,depth,time,var,varname,var_label,units,model,fpath,grid_type)
    % calculate weights
    rE = 6.371e6; % radius of Earth (m)
    dlat = lat(3) - lat(2); % spacing between latitudes
    area_3d = single(nan(length(lon),length(lat),length(depth)));
    depth_3d = repmat(permute(depth,[3 2 1]),length(lon),length(lat),1);
    for i = 1:length(lat)
        lat_area = 2*pi*rE^2*abs(sind(lat(i)-dlat/2)-sind(lat(i)+dlat/2));
        area_3d(:,i,:) = lat_area/length(lon); % m^2
    end
    vol = area_3d.*depth_3d; % m^3
    % calculate gloabl mean
    var_mean = nan(length(time),1);
    for t = 1:length(time)
        var_tmp = var(:,:,:,t);
        idx = ~isnan(var_tmp);
        var_mean(t) = sum(var_tmp(idx).*vol(idx))./sum(vol(idx));
    end
    % plot global mean
    figure('visible','off');
    set(gcf,'position',[100 100 800 400]);
    title(model);
    plot(double(time),var_mean,'LineWidth',3);
    datetick('x','keeplimits');
    ylabel([var_label ' (' units ')']);
    % save global mean plot
    if ~isfolder(['Figures/' model]); mkdir(['Figures/' model]); end
    export_fig(['Figures/' model '/' varname '_global_mean.png'],'-transparent');
    close
end
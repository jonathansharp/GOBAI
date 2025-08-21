function process_cmip_model(model,fpath,snap_date,start_year,rlz,grid_label,grid_type)

%% process date
start_month = 1;
date_str = num2str(snap_date);
end_year = str2double(date_str(1:4));
end_month = str2double(date_str(5:6));

%% define paths
path1_hist = [fpath 'historical/' grid_type '/'];
path1_ssp = [fpath 'ssp245/' grid_type '/'];
path2 = ['_Omon_' model '_'];
path3 = ['_' rlz '_' grid_label];
if ~isfolder([fpath 'combined/regridded']); mkdir([fpath 'combined/regridded']); end
% time extensions for cmip model paths
if strcmp(model,'GFDL-ESM4') || strcmp(model,'MPI-ESM1-2-LR')
    ext_hist1 = '199001-200912';
    ext_hist2 = '201001-201412';
    ext_hist3 = '';
    ext_ssp1 = '201501-203412';
    ext_ssp2 = ''; 
elseif strcmp(model,'NorESM-LM')
    ext_hist1 = '199001-199912';
    ext_hist2 = '200001-200912';
    ext_hist3 = '201001-201412';
    ext_ssp1 = '201501-202012';
    ext_ssp2 = '202101-203012';    
elseif strcmp(model,'CanESM5')
    ext_hist1 = '199101-200012';
    ext_hist2 = '200101-201012';
    ext_hist3 = '201101-201412';
    ext_ssp1 = '201501-202012';
    ext_ssp2 = '202101-203012';    
elseif strcmp(model,'ACCESS-ESM1-5')
    ext_hist1 = '199001-199912';
    ext_hist2 = '200001-200912';
    ext_hist3 = '201001-201412';
    ext_ssp1 = '201501-202412';
    ext_ssp2 = '202501-203412';
elseif strcmp(model,'IPSL-CM6A-LR')
    ext_hist1 = '195001-201412';
    ext_hist2 = '';
    ext_hist3 = '';
    ext_ssp1 = '201501-210012';
    ext_ssp2 = '';
end

%% load and process historical time
time_inf_historical = ncinfo([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'time');
% define origin date
origin_idx = find(strcmp({time_inf_historical.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_historical.Attributes(origin_idx).Value,'days since ')));
% read historical time
time_hist1 = ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'time');
if strcmp(model,'IPSL-CM6A-LR')
    time_hist2 = []; else
time_hist2 = ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist2 '.nc'],'time');
end
if strcmp(model,'GFDL-ESM4') || strcmp(model,'IPSL-CM6A-LR') || strcmp(model,'MPI-ESM1-2-LR')
    time_hist3 = []; else
time_hist3 = ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist3 '.nc'],'time');
end
% convert to matlab datenum
calendar_idx = find(strcmp({time_inf_historical.Attributes.Name},'calendar'));
if strcmp(time_inf_historical.Attributes(calendar_idx).Value,'365_day') || ...
    strcmp(time_inf_historical.Attributes(calendar_idx).Value,'noleap')
    time_hist = daynoleap2datenum([time_hist1;time_hist2;time_hist3]-1,origin_date(1));
elseif strcmp(time_inf_historical.Attributes(calendar_idx).Value,'gregorian') || ...
    strcmp(time_inf_historical.Attributes(calendar_idx).Value,'proleptic_gregorian')
    time_hist = datenum(origin_date)+datenum([time_hist1;time_hist2;time_hist3]-1);
end

%% load and process ssp time
time_inf_ssp = ncinfo([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp1 '.nc'],'time');
% define origin date
origin_idx = find(strcmp({time_inf_ssp.Attributes.Name},'units'));
origin_date = datevec(datenum(extractAfter(time_inf_ssp.Attributes(origin_idx).Value,'days since ')));
% read ssp time
time_ssp1 = ncread([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp1 '.nc'],'time');
if strcmp(model,'NorESM2-LM') || strcmp(model,'CanESM5') || strcmp(model,'ACCESS-ESM1-5') 
    time_ssp2 = ncread([path1_ssp 'o2' path2 'ssp245' path3 '_' ext_ssp2 '.nc'],'time');
else; time_ssp2 = []; end
% convert to matlab datenum
calendar_idx = find(strcmp({time_inf_ssp.Attributes.Name},'calendar'));
if strcmp(time_inf_ssp.Attributes(calendar_idx).Value,'365_day') || ...
    strcmp(time_inf_ssp.Attributes(calendar_idx).Value,'noleap')
    time_ssp = daynoleap2datenum([time_ssp1;time_ssp2]-1,origin_date(1));
elseif strcmp(time_inf_ssp.Attributes(calendar_idx).Value,'gregorian') || ...
    strcmp(time_inf_ssp.Attributes(calendar_idx).Value,'proleptic_gregorian')
    time_ssp = datenum(origin_date)+datenum([time_ssp1;time_ssp2]-1);
end

%% create time indices
time_min = datenum([start_year start_month 15]);
time_strt = find(abs(time_hist-time_min)==min(abs(time_hist-time_min)));
time_max = datenum([end_year end_month 15]);
time_end = find(abs(time_ssp-time_max)==min(abs(time_ssp-time_max)));
time = single([time_hist(time_strt:end);time_ssp(1:time_end)]);

%% load other 1-D model variables
if strcmp(grid_label,'gr')
    if strcmp(model,'GFDL-ESM4')
        lon = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lon'));
        lat = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lat'));
        lon2d = nan; lat2d = nan;
        depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lev'));
    end
elseif strcmp(grid_label,'gn')
    if strcmp(model,'NorESM2-LM') || strcmp(model,'IPSL-CM6A-LR')
        lon = (0.5:359.5)'; lat = (-89.5:89.5)';
        lon2d = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'nav_lon'));
        lon2d = convert_lon(lon2d,'0-360');
        lat2d = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'nav_lat'));
        depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'olevel'));
    elseif strcmp(model,'CanESM5') || strcmp(model,'ACCESS-ESM1-5') || strcmp(model,'MPI-ESM1-2-LR')
        lon = (0.5:359.5)'; lat = (-89.5:89.5)';
        lon2d = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'longitude'));
        lat2d = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'latitude'));
        depth = single(ncread([path1_hist 'o2' path2 'historical' path3 '_' ext_hist1 '.nc'],'lev'));
    end
end

% index to above 2100m
idx_depth = find(depth < 2100); depth = depth(idx_depth);

%% load, combine, and save potential temperature
nc_filepath = [fpath 'combined/regridded/thetao' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
nc_filepath_temp = [fpath 'combined/regridded/thetao' path2 ... % define temporary filepath
    'combined_' rlz '_gr-'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,nc_filepath_temp,path1_hist,path1_ssp,path2,path3,lon,lat,lon2d,lat2d,...
        depth,time,time_strt,grid_label,'thetao',idx_depth,'Sea Water Potential Temperature',...
        'degC',ext_hist1,ext_hist2,ext_hist3,ext_ssp1,ext_ssp2);
    % plot global mean timeseries
    thetao = ncread(nc_filepath,'thetao');
    plot_global_timeseries(lat,lon,depth,time,thetao,'thetao',...
        'Potential Temperature',[char(176) 'C'],model);
    clear thetao
end

%% load, combine, and save salinity
nc_filepath = [fpath 'combined/regridded/so' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
nc_filepath_temp = [fpath 'combined/regridded/so' path2 ... % define temporary filepath
    'combined_' rlz '_gr-'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,nc_filepath_temp,path1_hist,path1_ssp,path2,path3,lon,lat,lon2d,lat2d,...
        depth,time,time_strt,grid_label,'so',idx_depth,'Sea Water Salinity','n/a',...
        ext_hist1,ext_hist2,ext_hist3,ext_ssp1,ext_ssp2);
    % plot global mean timeseries
    so = ncread(nc_filepath,'so');
    plot_global_timeseries(lat,lon,depth,time,so,'so','Salinity','n/a',model);
    clear so
end

%% load, combine, and save chlorophyll
% nc_filepath = [fpath 'combined/regridded/chl' path2 ... % define filepath
%     'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
% nc_filepath_temp = [fpath 'combined/regridded/chl' path2 ... % define temporary filepath
%     'combined_' rlz '_gr-'];
% if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% if l ~= length(time)
%     combine_cmip_field(nc_filepath,nc_filepath_temp,path1_hist,path1_ssp,lon,lat,...
%         depth,time,time_strt,grid_label,'chl',idx_depth);
%     limit chlorophyll to only surface values (to match observations)
%     chl = repmat(chl(:,:,1,:),1,1,length(depth),1);
%     % plot global mean timeseries
%     chl = ncread(nc_filepath,'chl');
%     plot_global_timeseries(lat,lon,depth,time,chl,'chl','Chlorophyll','mg/m2',model,fpath,grid_type);
% end

%% process dimensional variables

%%%%%%%%%%%%%%%%%%%%%%%% GOTTA DO SOMETHING ABOUT 2D and 1D LONGITUDE

% establish dimensions
xdim = length(lon); ydim = length(lat); zdim = length(depth);
% expand coordinates to 3-D variables
lon_3d = repmat(lon,1,ydim,zdim);
lat_3d = repmat(lat',xdim,1,zdim);
depth_3d = repmat(permute(depth,[3 2 1]),xdim,ydim,1);
% calculate 3d pressure
if strcmp(model,'MPI-ESM1-2-LR')
    pres_3d = gsw_p_from_z(-depth_3d,lat_3d);
else
    pres_3d = -gsw_p_from_z(depth_3d,lat_3d);
end

%% calculate and save absolute salinity
nc_filepath = [fpath 'combined/regridded/abs_sal' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    so = ncread([fpath 'combined/regridded/so' path2 'combined_' ...
        rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'so');
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
        'Absolute Salinity','n/a',model);
    clear abs_sal
end

%% calculate and save conservative temperature
nc_filepath = [fpath 'combined/regridded/cns_tmp' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/regridded/thetao' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/regridded/abs_sal' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
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
        'Conservative Temperature',[char(176) 'C'],model);
    clear cns_tmp
end

%% calculate and save in situ temperature
nc_filepath = [fpath 'combined/regridded/tmp' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    thetao = ncread([fpath 'combined/regridded/thetao' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'thetao');
    abs_sal = ncread([fpath 'combined/regridded/abs_sal' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
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
        'In Situ Temperature',[char(176) 'C'],model);
    clear tmp
end

%% calculate and save potential density anomaly
nc_filepath = [fpath 'combined/regridded/sigma' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/regridded/cns_tmp' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/regridded/abs_sal' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
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
        'Potential Density','kg/m^3',model);
    clear sigma
end

%% calculate and save in situ density
nc_filepath = [fpath 'combined/regridded/dens' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    cns_tmp = ncread([fpath 'combined/regridded/cns_tmp' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'cns_tmp');
    abs_sal = ncread([fpath 'combined/regridded/abs_sal' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'abs_sal');
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
        'Density','kg/m^3',model);
    clear dens
end

%% calculate and save oxygen saturation
nc_filepath = [fpath 'combined/regridded/o2_sat' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
if l ~= length(time)
    tmp = ncread([fpath 'combined/regridded/tmp' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'tmp');
    so = ncread([fpath 'combined/regridded/so' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'so');
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
        'Oxygen Saturation','umol/kg',model);
    clear o2_sat
end

%% load, combine, and save o2
nc_filepath = [fpath 'combined/regridded/o2' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
nc_filepath_temp = [fpath 'combined/regridded/o2' path2 ... % define temporary filepath
    'combined_' rlz '_gr-'];
% check if file exists and is the proper size
if isfile(nc_filepath); l = nc_dim_length(nc_filepath,'time'); else; l = 0; end
% load and combine cmip output fields
if l ~= length(time)
    combine_cmip_field(model,nc_filepath,nc_filepath_temp,path1_hist,path1_ssp,path2,path3,lon,lat,lon2d,lat2d,...
        depth,time,time_strt,grid_label,'o2',idx_depth,'Dissolved Oxygen Content',...
        'mL/L',ext_hist1,ext_hist2,ext_hist3,ext_ssp1,ext_ssp2);
    o2 = ncread(nc_filepath,'o2');
    dens = ncread([fpath 'combined/regridded/dens' path2 ...
        'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'],'dens');
    o2 = (o2./dens).*(10^6); % convert oxygen from mol/m^3 to umol/kg
    ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
        {'lat' lat 'latitude' 'degrees north'},...
        {'depth' depth 'depth' 'meters'},...
        {'time' time 'time' 'days since 0000-01-01'},...
        {'o2' o2 'Dissolved Oxygen Content','umol/kg'});
    % plot global mean timeseries
    o2 = ncread(nc_filepath,'o2');
    plot_global_timeseries(lat,lon,depth,time,o2,'o2',...
        'Oxygen Amount Content','umol/kg',model);
    clear o2
end

%% add 3d pressure to NetCDFs
vars = {'o2' 'thetao' 'so' 'abs_sal' 'cns_tmp' 'tmp' 'sigma' 'dens'};
for v = 1:length(vars)
    if isfile([fpath 'combined/regridded/' vars{v} path2 'combined_' rlz '_regridded.nc'])
        if ~nc_var_exist([fpath 'combined/regridded/' vars{v} path2 'combined_' rlz '_regridded.nc'],'pres')
            nccreate([fpath 'combined/regridded/' vars{v} path2 'combined_' rlz '_regridded.nc'],...
                'pres','Dimensions',{'lon' length(lon) 'lat' length(lat) 'depth' length(depth)});
            ncwrite([fpath 'combined/regridded/' vars{v} path2 'combined_' rlz '_regridded.nc'],...
                'pres',pres_3d);
        end
    end
end

%% Plot T at 20 m
nc_filepath = [fpath 'combined/regridded/tmp' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(min(abs(depth-20))==abs(depth-20));
tmp = ncread(nc_filepath,'tmp',[1 1 idx_depth(1) 1],[Inf Inf 1 Inf]);
mean_tmp = mean(tmp,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
title(['Annual mean at ' num2str(depth(idx_depth(1))) ' m (' model ')'],'fontsize',16)
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
nc_filepath = [fpath 'combined/regridded/so' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(min(abs(depth-20))==abs(depth-20));
so = ncread(nc_filepath,'so',[1 1 idx_depth(1) 1],[Inf Inf 1 Inf]);
mean_so = mean(so,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
title(['Annual mean at ' num2str(depth(idx_depth(1))) ' m (' model ')'],'fontsize',16)
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
nc_filepath = [fpath 'combined/regridded/o2' path2 ... % define filepath
    'combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
idx_depth = find(min(abs(depth-20))==abs(depth-20));
o2 = ncread(nc_filepath,'o2',[1 1 idx_depth(1) 1],[Inf Inf 1 Inf]);
mean_o2 = mean(o2,4,'omitnan');
figure('visible','off');
worldmap([-90 90],[20 380]);
title(['Annual mean at ' num2str(depth(idx_depth(1))) ' m (' model ')'],'fontsize',16)
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
function combine_cmip_field(model,nc_filepath,nc_filepath_temp,path1_hist,path1_ssp,path2,path3,...
    lon,lat,lon2d,lat2d,depth,time,time_strt,grid_label,var_name,idx_depth,long_name,units,...
    ext_hist1,ext_hist2,ext_hist3,ext_ssp1,ext_ssp2)

if strcmp(model,'GFDL-ESM4') || strcmp(model,'MPI-ESM1-2-LR')

    % for models with two historical files and one ssp file
    dims_hist1 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],var_name);
    dims_hist2 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist2 '.nc'],var_name);

    % set up parallel pool with 20 workers
    tic; parpool; fprintf('Pool initiation: '); toc;

    parfor t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist1.Size(4)-time_strt
            idx_t = (time_strt-1)+t;
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)
            idx_t = time_strt+(t-dims_hist1.Size(4));
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist2 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist2.Size(1) dims_hist2.Size(2) max(idx_depth) 1]);
        else
            idx_t = time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4));
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        end
        % save variable to combined file
        save_temp_files(nc_filepath_temp,t,var,var_name,lon2d,lat2d,...
            lon,lat,depth,time,grid_label,long_name,units);
    end

elseif strcmp(model,'NorESM2-LM') || strcmp(model,'ACCESS-ESM1-5') || strcmp(model,'CanESM5')

    % for models with three historical files and two ssp files

    dims_hist1 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],var_name);
    dims_hist2 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist2 '.nc'],var_name);
    dims_hist3 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist3 '.nc'],var_name);
    dims_ssp1 = ncinfo([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp1 '.nc'],var_name);
    dims_ssp2 = ncinfo([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp2 '.nc'],var_name);

    % set up parallel pool with 20 workers
    tic; parpool; fprintf('Pool initiation: '); toc;

    parfor t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist1.Size(4)-time_strt
            idx_t = (time_strt-1)+t;
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)
            idx_t = time_strt+(t-dims_hist1.Size(4));
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist2 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist2.Size(1) dims_hist2.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)+dims_hist3.Size(4)
            idx_t = time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4));
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist3 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist3.Size(1) dims_hist3.Size(2) max(idx_depth) 1]);
        elseif t <= dims_hist1.Size(4)-time_strt+dims_hist2.Size(4)+dims_hist3.Size(4)+dims_ssp1.Size(4)
            idx_t = time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4)-dims_hist3.Size(4));
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_ssp1.Size(1) dims_ssp1.Size(2) max(idx_depth) 1]);
        else
            idx_t = time_strt+(t-dims_hist1.Size(4)-dims_hist2.Size(4)-dims_hist3.Size(4)-dims_ssp1.Size(4));
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp2 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_ssp2.Size(1) dims_ssp2.Size(2) max(idx_depth) 1]);
        end
        % save variable to combined file
        save_temp_files(nc_filepath_temp,t,var,var_name,lon2d,lat2d,...
            lon,lat,depth,time,grid_label,long_name,units);
    end

elseif strcmp(model,'IPSL-CM6A-LR')

    % for models with one historical file and one ssp file
    dims_hist1 = ncinfo([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],var_name);

    % set up parallel pool with 20 workers
    tic; parpool; fprintf('Pool initiation: '); toc;

    parfor t = 1:length(time)
        % read historical or ssp variable
        if t <= dims_hist1.Size(4)-time_strt
            idx_t = (time_strt-1)+t;
            var = ncread([path1_hist var_name path2 'historical' path3 '_' ext_hist1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
        else
            idx_t = time_strt+(t-dims_hist1.Size(4));
            var = ncread([path1_ssp var_name path2 'ssp245' path3 '_' ext_ssp1 '.nc'],...
                var_name,[1 1 1 idx_t],[dims_hist1.Size(1) dims_hist1.Size(2) max(idx_depth) 1]);
         end
        % save variable to combined file
        save_temp_files(nc_filepath_temp,t,var,var_name,lon2d,lat2d,...
            lon,lat,depth,time,grid_label,long_name,units);
    end

end

% end parallel session
delete(gcp('nocreate'));

% combine temporary files in main file
files = dir([nc_filepath_temp '*.nc']); % count files in folder
for t = 1:length(files)
    filename_temp = [nc_filepath_temp num2str(t) '.nc'];
    var_temp = ncread(filename_temp,var_name); % read
    if t == 1
        ncsave_4d(nc_filepath,{'lon' lon 'longitude' 'degrees east'},...
            {'lat' lat 'latitude' 'degrees north'},...
            {'depth' depth 'depth' 'meters'},...
            {'time' time(t) 'time' 'days since 0000-01-01'},...
            {var_name var_temp long_name units});
    else
        ncwrite(nc_filepath,'time',time(t),t);
        ncwrite(nc_filepath,var_name,var_temp,[1 1 1 t]);
    end
    % delete temporary file
    delete(filename_temp);
end

end

%% save variable to combined file
function save_temp_files(nc_filepath_temp,t,var,var_name,lon2d,lat2d,...
    lon,lat,depth,time,grid_label,long_name,units)

    % interpolate
    if strcmp(grid_label,'gr')
        var_interp = var;
    elseif strcmp(grid_label,'gn')
        [lon_interp,lat_interp] = ndgrid(lon,lat);
        var_interp = nan([size(lon_interp),length(depth)]);
        for d = 1:length(depth)
            warning('off','all'); % suppress warnings
            var_tmp = var(:,:,d); idx = ~isnan(var_tmp);
            interp_mod = scatteredInterpolant(double(lon2d(idx)),...
                double(lat2d(idx)),var_tmp(idx),'linear','linear');
            var_interp_tmp = interp_mod(lon_interp,lat_interp);
            interp_mod = scatteredInterpolant(double(lon2d(:)),...
                double(lat2d(:)),double(idx(:)),'nearest','nearest');
            mask_tmp = logical(interp_mod(lon_interp,lat_interp));
            var_interp_tmp(~mask_tmp) = NaN;
            var_interp(:,:,d) = var_interp_tmp;
            warning('on','all');
       end
    end

    % Write output in temporary files
    filename = [nc_filepath_temp num2str(t) '.nc'];
    if exist(filename,'file')==2; delete(filename); end
    nccreate(filename,'time','Dimensions',{'time' 1});
    ncwrite(filename,'time',time(t));
    nccreate(filename,var_name,'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'depth' length(depth)});
    ncwrite(filename,var_name,var_interp);


end

%% plot global timeseries
function plot_global_timeseries(lat,lon,depth,time,var,varname,var_label,units,model)
    % calculate volume
    vol = single(calculate_volume(lat,lon,depth));
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
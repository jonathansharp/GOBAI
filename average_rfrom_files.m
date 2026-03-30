%% average RFROM files into 1x1 degree monthly files
function average_rfrom_files(fpaths,ver,year,base_grid,start_year,end_year)

%% concatenate gobai in monthly files to match RFROM, and 1x1 degree files

% define gobai monthly mean file name
filename_monthly_mean_temp = [fpaths.temp_path 'RFROM_TEMP_' ver '/' ...
    'RFROMV' ver(2) ver(4) '_TEMP_MONTHLY_MEAN_1x1.nc'];
filename_monthly_mean_sal = [fpaths.sal_path 'RFROM_SAL_' ver '/' ...
    'RFROMV' ver(2) ver(4) '_SAL_MONTHLY_MEAN_1x1.nc'];
if isfile(filename_monthly_mean_temp); delete(filename_monthly_mean_temp); end
if isfile(filename_monthly_mean_sal); delete(filename_monthly_mean_sal); end
% new lat/lon
lon_new = (0.5:359.5)'; lon_bins_new = (0:365)';
lat_new = (-89.5:89.5)'; lat_bins_new = (-90:90)';
lon = ncread([fpaths.temp_path ...
    'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_1993_01.nc'],'longitude');
lat = ncread([fpaths.temp_path ...
    'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_1993_01.nc'],'latitude');
pres = ncread([fpaths.temp_path ...
    'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_1993_01.nc'],'mean_pressure');
prs_bnds = ncread([fpaths.temp_path ...
    'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_1993_01.nc'],'mean_pressure_bnds');

% read and write rfrom temperature schema
rfrom_info = ncinfo([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_1993_01.nc']);
% adjust variable-specific schema
rfrom_info.Variables(1).Size = length(lon_new);
rfrom_info.Variables(1).Dimensions.Length = length(lon_new);
rfrom_info.Variables(2).Size = length(lat_new);
rfrom_info.Variables(2).Dimensions.Length = length(lat_new);
rfrom_info.Variables(3).Size = Inf;
rfrom_info.Variables(3).Dimensions.Length = Inf;
for v = 1:length(rfrom_info.Variables)
    if strcmp(rfrom_info.Variables(v).Name,'ocean_temperature')
        var_idx = v;
    end
end
rfrom_info.Variables(var_idx).Size = [length(lon_new) length(lat_new) length(pres) Inf];
rfrom_info.Variables(var_idx).Dimensions(1).Length = length(lon_new);
rfrom_info.Variables(var_idx).Dimensions(2).Length = length(lat_new);
rfrom_info.Variables(var_idx).Dimensions(4).Length = Inf;
% adjust attributes
rfrom_info.Attributes(2).Value = ['RFROM ' ver ' (Monthly-Mean-1x1)'];
% adjust dimensions
rfrom_info.Dimensions(1).Length = length(lon_new);
rfrom_info.Dimensions(2).Length = length(lat_new);
rfrom_info.Dimensions(3).Length = Inf;
% write schema
ncwriteschema(filename_monthly_mean_temp,rfrom_info);

% read and write rfrom salinity schema
rfrom_info = ncinfo([fpaths.sal_path 'RFROM_SAL_v2.2/RFROMV22_SAL_STABLE_1993_01.nc']);
% adjust variable-specific schema
rfrom_info.Variables(1).Size = length(lon_new);
rfrom_info.Variables(1).Dimensions.Length = length(lon_new);
rfrom_info.Variables(2).Size = length(lat_new);
rfrom_info.Variables(2).Dimensions.Length = length(lat_new);
rfrom_info.Variables(3).Size = Inf;
rfrom_info.Variables(3).Dimensions.Length = Inf;
for v = 1:length(rfrom_info.Variables)
    if strcmp(rfrom_info.Variables(v).Name,'ocean_salinity')
        var_idx = v;
    end
end
rfrom_info.Variables(var_idx).Size = [length(lon_new) length(lat_new) length(pres) Inf];
rfrom_info.Variables(var_idx).Dimensions(1).Length = length(lon_new);
rfrom_info.Variables(var_idx).Dimensions(2).Length = length(lat_new);
rfrom_info.Variables(var_idx).Dimensions(4).Length = Inf;
% adjust attributes
rfrom_info.Attributes(2).Value = ['RFROM ' ver ' (Monthly-Mean-1x1)'];
% adjust dimensions
rfrom_info.Dimensions(1).Length = length(lon_new);
rfrom_info.Dimensions(2).Length = length(lat_new);
rfrom_info.Dimensions(3).Length = Inf;
% write schema
ncwriteschema(filename_monthly_mean_sal,rfrom_info);

% write dimensions to monthly mean file
ncwrite(filename_monthly_mean_temp,'longitude',lon_new);
ncwrite(filename_monthly_mean_temp,'latitude',lat_new);
ncwrite(filename_monthly_mean_temp,'mean_pressure',pres);
ncwrite(filename_monthly_mean_temp,'mean_pressure_bnds',prs_bnds);
ncwrite(filename_monthly_mean_sal,'longitude',lon_new);
ncwrite(filename_monthly_mean_sal,'latitude',lat_new);
ncwrite(filename_monthly_mean_sal,'mean_pressure',pres);
ncwrite(filename_monthly_mean_sal,'mean_pressure_bnds',prs_bnds);

% load rfrom dims
TS = load_RFROM_dim(fpaths.temp_path,'v2.2',start_year,end_year);

% loop through months
for m = 1:length(TS.months)
    % define temporary file names
    filename_temp = [fpaths.temp_path 'RFROM_TEMP_' ver '/RFROMV'...
        ver(2) ver(4) '_TEMP_STABLE_' num2str(TS.years(m)) '_' ...
        sprintf('%02d',TS.months(m)) '.nc'];
    filename_sal = [fpaths.sal_path 'RFROM_SAL_' ver '/RFROMV'...
        ver(2) ver(4) '_SAL_STABLE_' num2str(TS.years(m)) '_' ...
        sprintf('%02d',TS.months(m)) '.nc'];
    % read time
    time_temp = ncread(filename_temp,'time');
    time_sal = ncread(filename_sal,'time');
    % read temperature and salinity
    temp_4d = ncread(filename_temp,'ocean_temperature');
    sal_4d = ncread(filename_sal,'ocean_salinity');
    % write monthly mean 1x1 degree file
    temp_3d = mean(temp_4d,4,'omitnan');
    sal_3d = mean(sal_4d,4,'omitnan');  
    [lon_3d,lat_3d] = ndgrid(lon,lat);
    [~,~,Xnum] = histcounts(lon_3d(:),lon_bins_new);
    [~,~,Ynum] = histcounts(lat_3d(:),lat_bins_new);
    subs = [Xnum, Ynum];
    sz = [length(lon_new),length(lat_new)];
    temp_3d_1x1 = nan(length(lon_new),length(lat_new),length(pres));
    sal_3d_1x1 = nan(length(lon_new),length(lat_new),length(pres));
    for z = 1:length(pres)
        temp_3d_temp = temp_3d(:,:,z);
        temp_3d_1x1(:,:,z) = accumarray(subs,temp_3d_temp(:),sz,@nanmean);
        sal_3d_temp = sal_3d(:,:,z);
        sal_3d_1x1(:,:,z) = accumarray(subs,sal_3d_temp(:),sz,@nanmean);
    end
    % write temperature and salinity
    ncwrite(filename_monthly_mean_temp,'ocean_temperature',...
        temp_3d_1x1,[1 1 1 m]);
    ncwrite(filename_monthly_mean_sal,'ocean_salinity',...
        sal_3d_1x1,[1 1 1 m]);
    % write time to monthly mean 1x1 degree files
    ncwrite(filename_monthly_mean_temp,'time',mean(time_temp),m);
    ncwrite(filename_monthly_mean_sal,'time',mean(time_sal),m);
end
function subsample_cmip_model(param,param1,data,model,fpath,file_date,...
    snap_date,float_file_ext,start_year,grid_label,grid_type)

%% process date
date_str = num2str(snap_date);

% check if processed file already exists
if ~isfile([param1 '/Data/' model '_' param '_data_' file_date float_file_ext '.mat'])

%% define paths
path2 = ['_Omon_' model '_'];
path3 = ['_r1i1p1f1_' grid_label];

%% load variables
nc_filepath = [fpath 'combined/' grid_type '/o2' path2 ... % define filepath
    'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
lat = ncread(nc_filepath,'lat');
lon = ncread(nc_filepath,'lon');
time = ncread(nc_filepath,'time');
depth = ncread(nc_filepath,'depth');

%% bin data points
% put longitude on same scale (0 to 360)
data.longitude = convert_lon(data.longitude,'format','0:360');
lon = convert_lon(lon,'format','0:360');
% determine bin number of each longitude point
lon_diff = diff(lon)./2;
[~,~,Xnum] = histcounts(data.longitude,[lon(1)-lon_diff(1); ...
    lon(2:end)-lon_diff; lon(end)+lon_diff(end)]);
% determine bin number of each latitude point
lat_diff = diff(lat)./2;
[~,~,Ynum] = histcounts(data.latitude,[lat(1)-lat_diff(1); ...
    lat(2:end)-lat_diff; lat(end)+lat_diff(end)]);
% determine bin number of each latitude point
depth_diff = diff(depth)./2;
[~,~,Znum] = histcounts(data.depth,[depth(1)-depth_diff(1); ...
    depth(2:end)-depth_diff; depth(end)+depth_diff(end)]);
% determine bin number of each time point
time_diff = diff(time)./2;
[~,~,Tnum] = histcounts(data.time,[time(1)-time_diff(1); ...
    time(2:end)-time_diff; time(end)+time_diff(end)]);
clear lat_diff time_diff

%% accumulate 3D grid
idx = Xnum > 0 & Ynum > 0 & Znum > 0 & Tnum > 0;
subs = [Xnum(idx), Ynum(idx), Znum(idx) Tnum(idx)];
clear Xnum Ynum Znum 
sz = [length(lon),length(lat),length(depth),length(time)];
train_idx = accumarray(subs,data.oxygen(idx),sz);
% convert to train index
train_idx_mod = train_idx > 0;
clear train_idx

%% subsample model
vars = fieldnames(data);
for v = 1:length(vars)
    % define variable names from models
    if strcmp(vars{v},'oxygen'); var_name = 'o2';
    elseif strcmp(vars{v},'sigma'); var_name = 'sigma';
    elseif strcmp(vars{v},'temperature'); var_name = 'tmp';
    elseif strcmp(vars{v},'temperature_cns'); var_name = 'cns_tmp';
    elseif strcmp(vars{v},'salinity'); var_name = 'so';
    elseif strcmp(vars{v},'salinity_abs'); var_name = 'abs_sal';
    else var_name = 'skip';
    end
    if strcmp(var_name,'skip')
        % do nothing
    else
        nc_filepath = [fpath 'combined/' grid_type '/' var_name path2 ... % define filepath
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        var = ncread(nc_filepath,var_name);
        % zero trap for oxygen
        if strcmp(var_name,'o2'); var(var<0) = 0; end
        % extract model observations where we have data
        all_data.(vars{v}) = var(train_idx_mod);
    end
end

%% add 3d dimensional variables
lon = convert_lon(lon,'format','-180:180'); % convert longitude back to -180 to 180
% add longitude, latitude, and depth
[lon_3d,lat_3d,depth_3d] = ndgrid(lon,lat,depth);
lon_4d = repmat(lon_3d,1,1,1,length(time));
all_data.longitude = lon_4d(train_idx_mod); clear lon_4d
lat_4d = repmat(lat_3d,1,1,1,length(time));
all_data.latitude = lat_4d(train_idx_mod); clear lat_4d
depth_4d = repmat(depth_3d,1,1,1,length(time));
all_data.depth = depth_4d(train_idx_mod); clear depth_4d
% convert pressure from depth
all_data.pressure = gsw_p_from_z(-all_data.depth,all_data.latitude);

%% add time variables
time_4d = repmat(permute(time,[4 3 2 1]),length(lon),length(lat),length(depth),1);
all_data.time = time_4d(train_idx_mod); clear time_4d
% calculate day
date = datevec(double(all_data.time));
date0 = date;
date0(:,2:3) = 0;
all_data.day = datenum(date) - datenum(date0);
all_data.year = date(:,1);
% transform day by sine and cosine
all_data.day_sin = sin((2.*pi.*all_data.day)./365.25);
all_data.day_cos = cos((2.*pi.*all_data.day)./365.25);
% transform longitude by cosine:
all_data.lon_cos_1 = cosd(all_data.longitude-20);
all_data.lon_cos_2 = cosd(all_data.longitude-110);
% calculate bottom depth
% [lat_2d,lon_2d] = meshgrid(lat,lon);
% z = single(bottom_depth(lat_2d,lon_2d));

%% save data
save([param1 '/Data/' model '_' param '_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

else

%% display information
disp([model ' already subsampled for ' date_str(5:6) '/' date_str(1:4)]);

end

end
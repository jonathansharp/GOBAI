function subsample_cmip_model(param_props,data,model,fpath,file_date,...
    snap_date,float_file_ext,start_year,rlz)

%% process date
date_str = num2str(snap_date);

%% check if processed file already exists
if ~isfile([param_props.dir_name '/Data/' model '_' param_props.file_name '_data_' file_date float_file_ext '.mat'])
%% load variables
nc_filepath = [fpath 'combined/regridded/' param_props.file_name '_Omon_' model ... % define filepath
    '_combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
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
train_idx = accumarray(subs,data.(param_props.file_name)(idx),sz);
% convert to train index
train_idx_mod = train_idx > 0;
clear train_idx

%% subsample model
vars = fieldnames(data);
for v = 1:length(vars)
    % define variable names from models
    if strcmp(vars{v},'o2'); var_name = 'o2';
    elseif strcmp(vars{v},'no3'); var_name = 'no3';
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
        nc_filepath = [fpath 'combined/regridded/' var_name '_Omon_' model ... % define filepath
            '_combined_' rlz '_gr_' num2str(start_year) '01-' date_str '.nc'];
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

%% bin to lat/lon grid cells and plot gridded observations
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
% make plot
idx_depth = find(min(abs(all_data.depth-20))==abs(all_data.depth-20));
figure('visible','off'); hold on
worldmap([-90 90],[20 380]);
title('model');
title(['Annual mean at ' num2str(all_data.depth(idx_depth(1))) ' m (' model ')'],'fontsize',16)
set(gcf,'Position',[617, 599, 820, 420])
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
[lon_temp,z] = reformat_lon(lon,all_data.(['gridded_' param_props.file_name])(:,:,2),20);
pcolorm(lat,[lon_temp lon_temp(end)+1],[z;z(end,:)]');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([150 350]);
cmap = cmocean('ice'); cmap(1,:) = 1; colormap(cmap)
c.Label.String = ['Average Gridded ' param_props.fig_name];
c.FontSize = 12;
mlabel off; plabel off;
if ~isfolder(['Figures/' model]); mkdir(['Figures/' model]); end
export_fig(['Figures/' model '/subsampled_' param_props.fig_name '_20m.png'],'-transparent');
close

%% save data
save([param_props.dir_name '/Data/' model '_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
    'all_data','file_date','-v7.3');

else

%% display information
disp([model ' already subsampled for ' date_str(5:6) '/' date_str(1:4)]);

end

end
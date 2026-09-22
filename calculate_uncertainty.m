% calculate uncertainty
%
% DESCRIPTION:
% This function calculates uncertainty for GOBAI
% data products
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/28/2026

function calculate_uncertainty(param_props,base_grid,fpaths,...
    model_types,num_clusters,numWorkers_predict,file_date,snap_date,...
    float_file_ext,glodap_year,train_ratio,val_ratio,test_ratio,...
    flt,gld,osd,ctd)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; osd_ext = 'o'; else osd_ext = ''; end
if ctd == 1; ctd_ext = 'c'; else ctd_ext = ''; end

%% define gobai filepaths
gobai_alg_dir = ...
    [fpaths.param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
    '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
    '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/' ...
    float_ext glodap_ext osd_ext ctd_ext '/'];
gobai_filepath = [gobai_alg_dir 'GOBAI-' param_props.dir_name '-HR-v' ...
            num2str(snap_date) '.nc'];

%% load gobai dimensions
lon = ncread(gobai_filepath,'lon');
lat = ncread(gobai_filepath,'lat');
pres = ncread(gobai_filepath,'pres');
time = ncread(gobai_filepath,'time');
lon_0_360 = convert_lon(lon,'format','0-360');
[lon_0_360_3d,lat_3d,pres_3d] = ndgrid(lon_0_360,lat,pres);

%% calculate gridding uncertainty
% load data
load([param_props.dir_name '/Data/processed_all_' ...
    param_props.file_name '_data_' float_ext glodap_ext osd_ext ctd_ext ...
    '_' file_date float_file_ext '.mat'],'all_data');
% convert longitude
all_data.longitude = convert_lon(all_data.longitude,'format','0-360');
% calculate binned standard deviation
lon_diff = unique(diff(lon));
lon_edges = [lon - lon_diff;lon(end) + lon_diff];
lat_diff = unique(diff(lat));
lat_edges = [lat - lat_diff;lat(end) + lat_diff];
pres_edges = [0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000];
% get histogram counts in each bin
[~,~,Xnum] = histcounts(all_data.longitude,lon_edges);
[~,~,Ynum] = histcounts(all_data.latitude,lat_edges);
[~,~,Znum] = histcounts(all_data.pressure,pres_edges);
% calculate standard deviation in each bin
subs = [Xnum,Ynum,Znum];
sz = [length(lon),length(lat),length(pres)];
idx = ~any(subs==0,2);
stdev_grid = accumarray(subs,all_data.(param_props.file_name),sz,@nanstd);
sigma_grid = accumarray(subs,all_data.sigma,sz,@nanmean);
stdev_grid(stdev_grid==0) = NaN;
sigma_grid(sigma_grid==0) = NaN;
% calculate distance from coast
dist_3d = dist2coast(lat_3d,lon_0_360_3d);
% calculate bottom depth
bot_3d = bottom_depth(lat_3d,lon_0_360_3d);
% test plot
% figure; set(gcf,'position',[100 100 800 400]);
% pcolor(lon,lat,stdev_grid(:,:,11)'); shading flat; colorbar;
% caxis([0 30]); colormap(cmocean('balance','pivot',0));
% fit model of variability vs depth, sigma, and bottom depth
idx = ~isnan(stdev_grid) & ~isnan(pres_3d) & ~isnan(sigma_grid) & ~isnan(bot_3d);
y = log10(stdev_grid(idx));
X = [ones(size(y)) pres_3d(idx) dist_3d(idx) bot_3d(idx)];
[b,~,~,~,stats] = regress(y,X);
gobai_grid_uncer_main = 10.^(b(1) + b(2).*pres_3d + ...
    b(3).*dist_3d + b(4).*bot_3d);

%% calculate algorithm uncertainty
% loop through each model
for m = 1:length(model_types)
        
    % gobai filepath
    osse_filepath{m} = ...
        [fpaths.param_path 'GOBAI/' model_types{m} '/FFNN/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/train' ...
            num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
            num2str(100*test_ratio) '/' float_ext glodap_ext osd_ext ctd_ext ...
            '/GOBAI-' param_props.dir_name '-HR-v' num2str(snap_date) '.nc'];
    
    % delta filepath
    delta_filepath{m} = ...
        [fpaths.param_path 'GOBAI/' model_types{m} '/DELTA/c' ...
        num2str(num_clusters) '_' file_date float_file_ext '/' ...
        float_ext glodap_ext osd_ext ctd_ext];
    
    % load osse dimensions
    osse_lat = ncread(osse_filepath{m},'lat');
    osse_lon = ncread(osse_filepath{m},'lon');
    osse_depth = ncread(osse_filepath{m},'depth');
    osse_time = ncread(osse_filepath{m},'time');
    
    % calculate pressure
    [osse_lon_3d{m},osse_lat_3d{m},osse_depth_3d] = ...
        ndgrid(osse_lon,osse_lat,osse_depth);
    osse_pres_3d{m} = gsw_p_from_z(-abs(osse_depth_3d),osse_lat_3d{m});
    clear osse_depth_3d

end

% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

% loop through each timestep
parfor t = 1:length(time)
    
    % read in gobai for timestep and pre-define delta
    gobai = ncread(gobai_filepath,param_props.file_name,...
            [1 1 1 t],[Inf Inf Inf 1]);
    idx_gobai = ~isnan(gobai);
    gobai_delta = nan([size(gobai),length(model_types)]);
    
    % loop through each model
    for m = 1:length(model_types)
        
        % load difference between model and reconstruction
        delta_o2 = ncread([delta_filepath{m} '/delta_gobai-' ...
            param_props.file_name '.nc'],['delta_' param_props.file_name],...
            [1 1 1 t],[Inf Inf Inf 1]);
        
        % interpolate to gobai
        idx_del = ~isnan(delta_o2); % index to valid points for interpolation
        interp_temp = scatteredInterpolant(osse_lon_3d{m}(idx_del),...
            osse_lat_3d{m}(idx_del),osse_pres_3d{m}(idx_del),...
            delta_o2(idx_del),'linear','nearest'); % scattered interpolant with pressure
        gobai_delta_temp = nan(size(gobai));
        gobai_delta_temp(idx_gobai) = interp_temp(lon_0_360_3d(idx_gobai),...
            lat_3d(idx_gobai),pres_3d(idx_gobai));
        gobai_delta(:,:,:,m) = gobai_delta_temp;

    end

    % calculate algorithm uncertainty as root mean squared differences
    % between models and reconstructions
    gobai_alg_uncer = sqrt(mean(gobai_delta.^2,4,'omitnan'));

    % calculate measurement uncertainty as 3% of [O2]
    gobai_meas_uncer = gobai.*0.03;

    % calculate gridding uncertainty
    gobai_grid_uncer = gobai_grid_uncer_main;
    gobai_grid_uncer(~idx_gobai) = NaN;
    % figure; pcolor(lon,lat,gobai_grid_uncer(:,:,11)'); shading flat; colorbar; clim([0 15]);

    % save uncertainty for timestep in temporary files
    filename = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer-' num2str(t) '.nc'];
    if exist(filename,'file')==2; delete(filename); end
    % save time for timestep in temporary files
    nccreate(filename,'time','Dimensions',{'time' 1});
    ncwrite(filename,'time',time(t));
    % % algorithm uncertainty
    % nccreate(filename,['u_alg_' param_props.file_name],'Dimensions',...
    %     {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    % ncwrite(filename,['u_alg_' param_props.file_name],gobai_alg_uncer);
    % % measurment uncertainty
    % nccreate(filename,['u_meas_' param_props.file_name],'Dimensions',...
    %     {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    % ncwrite(filename,['u_meas_' param_props.file_name],gobai_meas_uncer);
    % % gridding uncertainty
    % nccreate(filename,['u_grid_' param_props.file_name],'Dimensions',...
    %     {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    % ncwrite(filename,['u_grid_' param_props.file_name],gobai_grid_uncer);
    % total uncertainty
    gobai_tot_uncer = sqrt(gobai_alg_uncer.^2 + gobai_meas_uncer.^2 + ...
        gobai_grid_uncer.^2);
    nccreate(filename,['u_tot_' param_props.file_name],'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    ncwrite(filename,['u_tot_' param_props.file_name],gobai_tot_uncer);

    % display information
    fprintf(['Algorithm Uncertainty Obtained for ' datestr(datenum(1950,1,1+time(t)))]);

end

% end parallel session
delete(gcp('nocreate'));

%% concatenate gobai uncertainty in main file
% create NetCDF that will be end result
filename = create_nc_file(gobai_alg_dir,param_props,lon,lat,pres,snap_date);
% loop through files
files = dir([gobai_alg_dir 'gobai-' param_props.file_name '-uncer-*.nc']); % count uncertainty files in folder
for cnt = 1:length(files)
    % define file name
    filename_temp = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer-' num2str(cnt) '.nc'];
    % read information from temporary file and write it to main file
    time = ncread(filename_temp,'time'); % read
    ncwrite(filename,'time',time,cnt); % write
    % gobai_3d = ncread(filename_temp,['u_alg_' param_props.file_name]); % read
    % ncwrite(filename,['u_alg_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    % gobai_3d = ncread(filename_temp,['u_meas_' param_props.file_name]); % read
    % ncwrite(filename,['u_meas_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    % gobai_3d = ncread(filename_temp,['u_grid_' param_props.file_name]); % read
    % ncwrite(filename,['u_grid_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    gobai_3d = ncread(filename_temp,['u_tot_' param_props.file_name]); % read
    ncwrite(filename,['u_tot_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    % delete temporary file
    delete(filename_temp);
end

end

%% for creating main netCDF file
function filename = create_nc_file(gobai_alg_dir,param_props,lon,lat,pres,snap_date)

xdim = length(lon);
ydim = length(lat);
zdim = length(pres);

% define file name
filename = [gobai_alg_dir 'GOBAI-' param_props.dir_name '-HR-v' ...
    num2str(snap_date) '-Uncer.nc'];

% create folder and file
if ~isfolder(gobai_alg_dir); mkdir(gobai_alg_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% % bgc parameter (alg)
% nccreate(filename,['u_alg_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
%     'DataType','single','FillValue',NaN);
% ncwriteatt(filename,['u_alg_' param_props.file_name],'units',param_props.units);
% % bgc parameter (meas)
% nccreate(filename,['u_meas_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
%     'DataType','single','FillValue',NaN);
% ncwriteatt(filename,['u_meas_' param_props.file_name],'units',param_props.units);
% % bgc parameter (grid)
% nccreate(filename,['u_grid_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
%     'DataType','single','FillValue',NaN);
% ncwriteatt(filename,['u_grid_' param_props.file_name],'units',param_props.units);
% bgc parameter (tot)
nccreate(filename,['u_tot_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,['u_tot_' param_props.file_name],'units',param_props.units);
% time
nccreate(filename,'time','Dimensions',{'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,'time','units','days since 1950-0-0');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');
% longitude
nccreate(filename,'lon','Dimensions',{'lon',xdim},...
    'DataType','single','FillValue',NaN);
ncwrite(filename,'lon',lon);
ncwriteatt(filename,'lon','units','degrees_east');
ncwriteatt(filename,'lon','axis','X');
ncwriteatt(filename,'lon','long_name','longitude');
ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
% latitude
nccreate(filename,'lat','Dimensions',{'lat',ydim},...
    'DataType','single','FillValue',NaN);
ncwrite(filename,'lat',lat);
ncwriteatt(filename,'lat','units','degrees_north');
ncwriteatt(filename,'lat','axis','Y');
ncwriteatt(filename,'lat','long_name','latitude');
ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
% pressure
nccreate(filename,'pres','Dimensions',{'pres',zdim},...
    'DataType','single','FillValue',NaN);
ncwrite(filename,'pres',pres);
ncwriteatt(filename,'pres','units','decibars');
ncwriteatt(filename,'pres','axis','Z');
ncwriteatt(filename,'pres','long_name','pressure');
ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');

end

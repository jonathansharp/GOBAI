% calculate uncertainty
%
% DESCRIPTION:
% This function calculates uncertainty for GOBAI
% data products
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/7/2024

function calculate_uncertainty(param_props,base_grid,param_path,...
    temp_path,sal_path,fpath,model_types,num_clusters,...
    numWorkers_predict,file_date,float_file_ext,glodap_year,...
    train_ratio,val_ratio,test_ratio)

%% define gobai filepaths
gobai_alg_dir = ...
    [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
    '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
    '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
gobai_filepath = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];

%% load gobai dimensions
lon = ncread(gobai_filepath,'lon');
lat = ncread(gobai_filepath,'lat');
pres = ncread(gobai_filepath,'pres');
time = ncread(gobai_filepath,'time');
lon_0_360 = convert_lon(lon,'format','0-360');
[lon_0_360_3d,lat_3d,pres_3d] = ndgrid(lon_0_360,lat,pres);

%% determine coefficients to calculate gridding uncertainty
[b,stats] = calculate_gridding_uncertainty_coeffs(file_date,...
    temp_path,sal_path,float_file_ext,glodap_year,lon,lat,time,lon_0_360_3d,lat_3d,pres_3d);

%% loop through each model
for m = 1:length(model_types)

        %% combined filepath
        % gobai
        osse_filepath{m} = ...
            [fpath model_types{m} '/GOBAI/' model_types{m} '/FFNN/c' ...
            num2str(num_clusters) '_' file_date float_file_ext '/train' ...
                num2str(100*train_ratio) '_val' num2str(100*val_ratio) '_test' ...
                num2str(100*test_ratio) '/gobai-' param_props.file_name '.nc'];
        % delta
        delta_filepath{m} = [fpath model_types{m} '/GOBAI/' model_types{m} '/DELTA/c' ...
            num2str(num_clusters) '_' file_date float_file_ext];
        
        %% load osse dimensions
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

%% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% loop through each month
parfor t = 1:length(time)

    % read in gobai for timestep and pre-define delta
    gobai = ncread(gobai_filepath,param_props.file_name,...
            [1 1 1 t],[Inf Inf Inf 1]);
    idx_gobai = ~isnan(gobai);
    gobai_delta = nan([size(gobai),length(model_types)]);

    %% loop through each model
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
    bot_3d = bottom_depth(lat_3d,lon_0_360_3d);
    bot_3d(bot_3d > 2000) = 2000;
    gobai_grid_uncer = b(1) + b(2).*pres_3d + b(3).*bot_3d;
    gobai_grid_uncer(~idx_gobai) = NaN;
    % figure; pcolor(lon,lat,gobai_grid_uncer(:,:,11)'); shading flat; colorbar; clim([0 15]);

    % save uncertainty for timestep in temporary files
    filename = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer-' num2str(t) '.nc'];
    if exist(filename,'file')==2; delete(filename); end
    % time
    nccreate(filename,'time','Dimensions',{'time' 1});
    ncwrite(filename,'time',time(t));
    % algorithm uncertainty
    nccreate(filename,['u_alg_' param_props.file_name],'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    ncwrite(filename,['u_alg_' param_props.file_name],gobai_alg_uncer);
    % measurment uncertainty
    nccreate(filename,['u_meas_' param_props.file_name],'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    ncwrite(filename,['u_meas_' param_props.file_name],gobai_meas_uncer);
    % gridding uncertainty
    nccreate(filename,['u_grid_' param_props.file_name],'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    ncwrite(filename,['u_grid_' param_props.file_name],gobai_grid_uncer);
    % total uncertainty
    gobai_tot_uncer = sqrt(gobai_alg_uncer.^2 + gobai_meas_uncer.^2 + ...
        gobai_grid_uncer.^2);
    nccreate(filename,['u_tot_' param_props.file_name],'Dimensions',...
        {'lon' length(lon) 'lat' length(lat) 'pres' length(pres)});
    ncwrite(filename,['u_tot_' param_props.file_name],gobai_tot_uncer);

    % display information
    % fprintf(['Algorithm Uncertainty Obtained for']);

end

% end parallel session
delete(gcp('nocreate'));

%% concatenate gobai uncertainty in main file
% create NetCDF that will be end result
create_nc_file(gobai_alg_dir,param_props,lon,lat,pres);
% 
files = dir([gobai_alg_dir 'gobai-' param_props.file_name '-uncer-*.nc']); % count uncertainty files in folder
filename = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer.nc'];
for cnt = 1:length(files)
    % define file name
    filename_temp = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer-' num2str(cnt) '.nc'];
    % read information from temporary file and write it to main file
    time = ncread(filename_temp,'time'); % read
    ncwrite(filename,'time',time,cnt); % write
    gobai_3d = ncread(filename_temp,['u_alg_' param_props.file_name]); % read
    ncwrite(filename,['u_alg_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    gobai_3d = ncread(filename_temp,['u_meas_' param_props.file_name]); % read
    ncwrite(filename,['u_meas_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    gobai_3d = ncread(filename_temp,['u_grid_' param_props.file_name]); % read
    ncwrite(filename,['u_grid_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    gobai_3d = ncread(filename_temp,['u_tot_' param_props.file_name]); % read
    ncwrite(filename,['u_tot_' param_props.file_name],gobai_3d,[1 1 1 cnt]); % write
    % delete temporary file
    delete(filename_temp);
end

end

%% subfunction for calculating coefficients for gridding uncertainty
function [b,stats] = calculate_gridding_uncertainty_coeffs(file_date,temp_path,sal_path,...
    float_file_ext,glodap_year,lon,lat,time,lon_0_360_3d,lat_3d,pres_3d)

%% load interpolated float and glodap data
load(['O2/Data/processed_float_o2_data_' file_date float_file_ext '.mat'],...
    'float_data');
load(['O2/Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% determine histogram counts and indices
% establish edges of bins
x_edges = min(floor(lon)):1:max(ceil(lon)); x_bins = lon;
y_edges = min(floor(lat)):1:max(ceil(lat)); y_bins = lat;
z_edges = [0 5:10:175 190:20:450 475:50:1375 1450:100:1950 2000];
z_bins = [2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975];
t_edges = [datenum(1949,12,16)+time(1);datenum(1950,1,15)+time];
t_bins = datenum(1950,0,0)+time;
% get histogram counts in each bin for each dataset
float_lon_conv = convert_lon(float_data.LON,'format','0-360');
float_lon_conv(float_lon_conv<20) = float_lon_conv(float_lon_conv<20) + 360;
glodap_lon_conv = convert_lon(glodap_data.LON,'format','0-360');
glodap_lon_conv(glodap_lon_conv<20) = glodap_lon_conv(glodap_lon_conv<20) + 360;
[~,~,Xnum_float] = histcounts(float_lon_conv,x_edges);
[~,~,Ynum_float] = histcounts(float_data.LAT,y_edges);
[~,~,Znum_float] = histcounts(float_data.PRES,z_edges);
[~,~,Tnum_float] = histcounts(float_data.TIME,t_edges);
[~,~,Xnum_glodap] = histcounts(glodap_lon_conv,x_edges);
[~,~,Ynum_glodap] = histcounts(glodap_data.LAT,y_edges);
[~,~,Znum_glodap] = histcounts(glodap_data.PRES,z_edges);
[~,~,Tnum_glodap] = histcounts(glodap_data.TIME,t_edges);
% accumulate index of counts
subs_float = [Xnum_float,Ynum_float,Znum_float,Tnum_float];
idx_float = ~any(subs_float==0,2);
subs_glodap = [Xnum_glodap,Ynum_glodap,Znum_glodap,Tnum_glodap];
idx_glodap = ~any(subs_glodap==0,2);
% determine size of 4D grid
sz = [length(x_bins),length(y_bins),length(z_bins),length(t_bins)];

% get variability within each grid cell (oxygen)
float_oxy_std = single(accumarray(subs_float(idx_float,:),float_data.OXY(idx_float),sz,@nanstd,nan));
glodap_oxy_std = single(accumarray(subs_glodap(idx_glodap,:),glodap_data.OXY(idx_glodap),sz,@nanstd,nan));
oxy_std = mean(cat(5,float_oxy_std,glodap_oxy_std),5,'omitnan');
clear float_oxy_std glodap_oxy_std
idx = ~isnan(oxy_std); % index when oxygen variability exists
oxy_std = double(oxy_std(idx));

% temperature and salinity variability
tmp = ncread([temp_path 'RG_Climatology_Temp.nc'],'Temperature');
tmp_var = repmat(double(std(tmp,[],4)),1,1,1,length(time)); clear tmp;
tmp_var_fit = tmp_var(idx);
sal = ncread([temp_path 'RG_Climatology_Sal.nc'],'Salinity');
sal_var = repmat(double(std(sal,[],4)),1,1,1,length(time)); clear sal;
sal_var_fit = sal_var(idx);

% get mean within each grid cell (sigma)
float_sig = single(accumarray(subs_float(idx_float,:),float_data.SIGMA(idx_float),sz,@nanmean,nan));
glodap_sig = single(accumarray(subs_glodap(idx_glodap,:),glodap_data.SIGMA(idx_glodap),sz,@nanmean,nan));
sig_fit = mean(cat(5,float_sig,glodap_sig),5,'omitnan');
clear float_sig glodap_sig
sig_fit = sig_fit(idx);

% calculate distance from coast
dist_fit = dist2coast(lat_3d,lon_0_360_3d);
dist_fit = repmat(dist_fit,1,1,1,length(time));
dist_fit = dist_fit(idx);

% replicate pressure
pres_fit = repmat(pres_3d,1,1,1,length(time));
pres_fit = pres_fit(idx);

% calculate bottom depth
bot_3d = bottom_depth(lat_3d,lon_0_360_3d);
bot_fit = repmat(bot_3d,1,1,1,length(time));
bot_fit = bot_fit(idx);

% remove standard deviations calculated from nine or fewer measurements
oxygen_count = accumarray(subs_float(idx_float,:),1,sz) + accumarray(subs_glodap(idx_glodap,:),1,sz);
oxygen_count = oxygen_count(idx);
oxygen_count_idx = oxygen_count < 10;
oxy_std(oxygen_count_idx) = [];
sig_fit(oxygen_count_idx) = [];
pres_fit(oxygen_count_idx) = [];
dist_fit(oxygen_count_idx) = [];
bot_fit(oxygen_count_idx) = [];
bot_fit(bot_fit>2000) = 2000;
tmp_var_fit(oxygen_count_idx) = [];
sal_var_fit(oxygen_count_idx) = [];
% scatter different parameters against std
% figure; scatter(sig_fit,oxy_std);
% fit model of variability vs depth and distance from shore
[b,~,~,~,stats] = ...
    regress(oxy_std,[ones(size(pres_fit)) pres_fit bot_fit]);
stats(1); % R^2 statistic

end

%% for creating main netCDF file
function create_nc_file(gobai_alg_dir,param_props,lon,lat,pres)

xdim = length(lon);
ydim = length(lat);
zdim = length(pres);

% define file name
filename = [gobai_alg_dir 'gobai-' param_props.file_name '-uncer.nc'];

% create folder and file
if ~isfolder(gobai_alg_dir); mkdir(gobai_alg_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% bgc parameter (alg)
nccreate(filename,['u_alg_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,['u_alg_' param_props.file_name],'units',param_props.units);
% bgc parameter (meas)
nccreate(filename,['u_meas_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,['u_meas_' param_props.file_name],'units',param_props.units);
% bgc parameter (grid)
nccreate(filename,['u_grid_' param_props.file_name],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,['u_grid_' param_props.file_name],'units',param_props.units);
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

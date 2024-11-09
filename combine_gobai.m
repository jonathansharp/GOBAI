% combine_gobai
%
% DESCRIPTION:
% This function combines output from gobai gridded
% fields obtained via different machine learning algorithms.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/22/2023

function combine_gobai(param,fpath,base_grid,file_date,float_file_ext,...
    num_clusters,start_year,snap_date,train_ratio,val_ratio,...
    test_ratio,numtrees,minLeafSize,numstumps,numbins,varargin)

%% set defaults and process optional input arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    elseif strcmpi(varargin{i}, 'grid_label')
        grid_label = varargin{i+1};
    elseif strcmpi(varargin{i}, 'grid_type')
        grid_type = varargin{i+1};
    end
end

%% process date
date_str = num2str(snap_date);

%% process parameter name
[param1,~,~,~,~,~,~,~,units,long_param_name] = param_name(param);

%% create directory names
gobai_ffnn_dir = ... % FFNN
    [param1 '/Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
gobai_rfr_dir = ... % RFR
    [param1 '/Data/GOBAI/' base_grid '/RFR/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];
gobai_gbm_dir = ... % GBM
    [param1 '/Data/GOBAI/' base_grid '/GBM/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numstumps) '_bin' num2str(numbins) '/'];
gobai_dir = ... % final product
    [param1 '/Data/GOBAI/' base_grid '/AVG/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/'];

%% create netCDF file
if strcmp(base_grid,'RG')
    TS = load_RG_dim([fpath '/Data/RG_CLIM/']);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param,...
        units,long_param_name);
elseif strcmp(base_grid,'RFROM')
    TS = load_RFROM_dim([fpath '/Data/RFROM/']);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param,...
        units,long_param_name);
else
    % define paths
    path2 = ['_Omon_' base_grid '_']; path3 = '_r1i1p1f1_gr';
    % define filepaths
    nc_filepath_abs_sal = [fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param,...
        units,long_param_name);
end

%% average over each time window
m_last = (str2num(date_str(1:4))-2004)*12+str2num(date_str(5:6));
for m = 1:12%m_last

    %% set counter
    cnt = m;

    %% load monthly outputs
    % load monthly output (FFNN)
    gobai_3d_ffnn = ncread([gobai_ffnn_dir 'gobai-' param '.nc'],...
        param,[1 1 1 cnt],[Inf Inf Inf 1]);
    % load monthly output (RFR)
    gobai_3d_rfr = ncread([gobai_rfr_dir 'gobai-' param '.nc'],...
        param,[1 1 1 cnt],[Inf Inf Inf 1]);
    % load monthly output (GBM)
    gobai_3d_gbm = ncread([gobai_gbm_dir 'gobai-' param '.nc'],...
        param,[1 1 1 cnt],[Inf Inf Inf 1]);

    %% average monthly outputs
    gobai_3d_avg = mean(cat(4,gobai_3d_ffnn,gobai_3d_rfr,gobai_3d_gbm),4,'omitnan');

    %% create folder and save monthly output
    filename = [gobai_dir 'gobai-' param '.nc'];
    % write to file
    ncwrite(filename,'time',datenum(TS.years(cnt),TS.months(cnt),15),cnt);
    ncwrite(filename,param,gobai_3d_avg,[1 1 1 cnt]);

end

end

%% function for creating netCDF file
function create_nc_file(TS,base_grid,xdim,ydim,zdim,gobai_dir,...
    param,units,long_param_name)

% define file name
filename = [gobai_dir 'gobai-' param '.nc'];

% create folder and file
if ~isfolder([pwd '/' gobai_dir]); mkdir(gobai_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% bgc parameter
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,param,'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
else
    nccreate(filename,param,'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
end
ncwriteatt(filename,param,'units',units);
ncwriteatt(filename,param,'long_name',long_param_name);
% longitude
nccreate(filename,'lon','Dimensions',{'lon',xdim},...
    'DataType','single','FillValue',NaN);
ncwrite(filename,'lon',TS.Longitude);
ncwriteatt(filename,'lon','units','degrees_east');
ncwriteatt(filename,'lon','axis','X');
ncwriteatt(filename,'lon','long_name','longitude');
ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
% latitude
nccreate(filename,'lat','Dimensions',{'lat',ydim},...
    'DataType','single','FillValue',NaN);
ncwrite(filename,'lat',TS.Latitude);
ncwriteatt(filename,'lat','units','degrees_north');
ncwriteatt(filename,'lat','axis','Y');
ncwriteatt(filename,'lat','long_name','latitude');
ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
% pressure (or depth)
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,'pres','Dimensions',{'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'pres',TS.Pressure);
    ncwriteatt(filename,'pres','units','decibars');
    ncwriteatt(filename,'pres','axis','Z');
    ncwriteatt(filename,'pres','long_name','pressure');
    ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');
else
    nccreate(filename,'depth','Dimensions',{'depth',zdim},...
    'DataType','single','FillValue',NaN);
    ncwrite(filename,'depth',TS.Depth);
    ncwriteatt(filename,'depth','units','meters');
    ncwriteatt(filename,'depth','axis','Z');
    ncwriteatt(filename,'depth','long_name','depth');
    ncwriteatt(filename,'depth','_CoordinateAxisType','Depth');
end
% time
nccreate(filename,'time','Dimensions',{'time',Inf},...
    'DataType','single','FillValue',NaN);
ncwriteatt(filename,'time','units','days since 0000-01-01');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');

end
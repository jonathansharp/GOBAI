% combine_gobai
%
% DESCRIPTION:
% This function combines output from gobai gridded
% fields obtained via different machine learning algorithms.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/10/2025

% function combine_gobai(param_props,temp_path,param_path,base_grid,file_date,float_file_ext,...
%     num_clusters,start_year,end_year,snap_date,train_ratio,val_ratio,...
%     test_ratio,numtrees,minLeafSize,numstumps,numbins,varargin)
function combine_gobai(param_props,temp_path,param_path,base_grid,file_date,float_file_ext,...
    num_clusters_1,num_clusters_2,num_clusters_3,start_year,end_year,snap_date,train_ratio,val_ratio,...
    test_ratio,numtrees,minLeafSize,numstumps,numbins,varargin)

%% set defaults and process optional input arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    end
end

%% process date
date_str = num2str(snap_date);

%% create directory names (multiple models)
% gobai_ffnn_dir = ... % FFNN
%     [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
%     float_file_ext '/train' num2str(100*train_ratio) '_val' ...
%     num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
% gobai_rfr_dir = ... % RFR
%     [param_path 'GOBAI/' base_grid '/RFR/c' num2str(num_clusters) '_' file_date ...
%     float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];
% gobai_gbm_dir = ... % GBM
%     [param_path 'GOBAI/' base_grid '/GBM/c' num2str(num_clusters) '_' file_date ...
%     float_file_ext '/tr' num2str(numstumps) '_bin' num2str(numbins) '/'];
% gobai_dir = ... % final product
%     [param_path 'GOBAI/' base_grid '/AVG/c' num2str(num_clusters) '_' file_date ...
%     float_file_ext '/'];

%% create directory names (multiple cluster formations)
gobai_ffnn_dir_1 = ... % FFNN
    [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters_1) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
gobai_ffnn_dir_2 = ... % FFNN
    [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters_2) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
gobai_ffnn_dir_3 = ... % FFNN
    [param_path 'GOBAI/' base_grid '/FFNN/c' num2str(num_clusters_3) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
gobai_dir = ... % final product
    [param_path 'GOBAI/' base_grid '/AVG/multiple_clusters_' file_date ...
    float_file_ext '/'];

%% create netCDF file
if strcmp(base_grid,'RG')
    TS = load_RG_dim(temp_path);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param_props);
elseif strcmp(base_grid,'RFROM')
    TS = load_RFROM_dim(temp_path,start_year,end_year);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param_props);
else
    % define paths
    path2 = ['_Omon_' base_grid '_']; path3 = '_r1i1p1f1_gr';
    % define filepaths
    nc_filepath_abs_sal = [temp_path 'combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_dir,param_props);
end

%% average over each time window
for m = 1:length(TS.Time)

        %% set counter
        cnt = m;

        % %% load monthly outputs (multiple models)
        % % load monthly output (FFNN)
        % gobai_3d_ffnn = ncread([gobai_ffnn_dir 'gobai-' param_props.file_name '.nc'],...
        %     param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
        % % load monthly output (RFR)
        % gobai_3d_rfr = ncread([gobai_rfr_dir 'gobai-' param_props.file_name '.nc'],...
        %     param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
        % % load monthly output (GBM)
        % gobai_3d_gbm = ncread([gobai_gbm_dir 'gobai-' param_props.file_name '.nc'],...
        %     param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
 
        %% load monthly outputs (multiple cluster formations)
        % load monthly output (FFNN)
        gobai_3d_ffnn_1 = ncread([gobai_ffnn_dir_1 'gobai-' param_props.file_name '.nc'],...
            param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
        gobai_3d_ffnn_2 = ncread([gobai_ffnn_dir_2 'gobai-' param_props.file_name '.nc'],...
            param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
        gobai_3d_ffnn_3 = ncread([gobai_ffnn_dir_3 'gobai-' param_props.file_name '.nc'],...
            param_props.file_name,[1 1 1 cnt],[Inf Inf Inf 1]);
        
        % %% average monthly outputs (multiple models)
        % gobai_3d_avg = mean(cat(4,gobai_3d_ffnn,gobai_3d_rfr,gobai_3d_gbm),4,'omitnan');
        % gobai_3d_var = std(cat(4,gobai_3d_ffnn,gobai_3d_rfr,gobai_3d_gbm),[],4,'omitnan');

        %% average monthly outputs (multiple cluster formations)
        gobai_3d_avg = mean(cat(4,gobai_3d_ffnn_1,gobai_3d_ffnn_2,gobai_3d_ffnn_3),4,'omitnan');
        gobai_3d_var = std(cat(4,gobai_3d_ffnn_1,gobai_3d_ffnn_2,gobai_3d_ffnn_3),[],4,'omitnan');
    
        %% create folder and save monthly output
        filename = [gobai_dir 'gobai-' param_props.file_name '.nc'];
        % write to file
        if strcmp(base_grid,'RFROM')
            ncwrite(filename,'time',TS.Time(m),cnt);
        else
            ncwrite(filename,'time',datenum(TS.years(cnt),TS.months(cnt),15)-datenum(1950,0,0),cnt);
        end
        ncwrite(filename,param_props.file_name,gobai_3d_avg,[1 1 1 cnt]);
        ncwrite(filename,[param_props.file_name '_var'],gobai_3d_var,[1 1 1 cnt]);

end

end

%% function for creating netCDF file
function create_nc_file(TS,base_grid,xdim,ydim,zdim,gobai_dir,param_props)

% define file name
filename = [gobai_dir 'gobai-' param_props.file_name '.nc'];

% create folder and file
if ~isfolder(gobai_dir); mkdir(gobai_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% bgc parameter
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
    nccreate(filename,[param_props.file_name '_var'],'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
else
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
    nccreate(filename,[param_props.file_name '_var'],'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
end
ncwriteatt(filename,param_props.file_name,'units',param_props.units);
ncwriteatt(filename,param_props.file_name,'long_name',param_props.long_param_name);
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
ncwriteatt(filename,'time','units','days since 1950-0-0');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');

end
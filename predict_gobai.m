% predict_gobai
%
% DESCRIPTION:
% This function uses trained models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/6/2024

function predict_gobai(alg_type,param,fpath,base_grid,file_date,float_file_ext,...
    num_clusters,variables,thresh,numWorkers_predict,start_year,...
    snap_date,varargin)

%% set defaults and process optional input arguments
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    elseif strcmpi(varargin{i}, 'grid_label')
        grid_label = varargin{i+1};
    elseif strcmpi(varargin{i}, 'grid_type')
        grid_type = varargin{i+1};
    elseif strcmpi(varargin{i}, 'train_ratio')
            train_ratio = varargin{i+1};
    elseif strcmpi(varargin{i}, 'val_ratio')
        val_ratio = varargin{i+1};
    elseif strcmpi(varargin{i}, 'test_ratio')
        test_ratio = varargin{i+1};
    elseif strcmpi(varargin{i}, 'numtrees')
        numtrees = varargin{i+1};
    elseif strcmpi(varargin{i}, 'minLeafSize')
        minLeafSize = varargin{i+1};
    elseif strcmpi(varargin{i}, 'numstumps')
        numstumps = varargin{i+1};
    elseif strcmpi(varargin{i}, 'numbins')
        numbins = varargin{i+1};
    end
end
% make sure alg_type is valid
if ~strcmp(alg_type,'FFNN') | ~strcmp(alg_type,'RFR') | ~strcmp(alg_type,'GBM')
    disp('"alg_type" must be "FFNN", "RFR", or "GBM"')
end

%% process date
date_str = num2str(snap_date);

%% directory base
if strcmp(alg_type,'FFNN')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
elseif strcmp(alg_type,'RFR')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
elseif strcmp(alg_type,'GBM')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps;numbins});
end

%% process parameter name
[param1,param2,~,~,~,~,~,~,units,long_param_name] = param_name(param);

%% create directory and file names
alg_dir = [param1 '/Models/' dir_base];
alg_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    alg_fnames(c) = ...
        {[alg_type '_' param2 '_C' num2str(c)]};
end
if strcmp(alg_type,'FFNN')
    gobai_alg_dir = ...
        [param1 '/Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
        '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];
elseif strcmp(alg_type,'RFR')
    gobai_alg_dir = ...
        [param1 '/Data/GOBAI/' base_grid '/RFR/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numtrees) '_lf' ...
        num2str(minLeafSize) '/'];
elseif strcmp(alg_type,'GBM')
    gobai_alg_dir = ...
        [param1 '/Data/GOBAI/' base_grid '/GBM/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numstumps) ...
        '_bin' num2str(numbins) '/'];
end

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% create netCDF file
if strcmp(base_grid,'RG')
    TS = load_RG_dim([fpath '/Data/RG_CLIM/']);
   % create file
   create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param,...
        units,long_param_name);
elseif strcmp(base_grid,'RFROM')
    TS = load_RFROM_dim([fpath '/Data/RFROM/']);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param,...
        units,long_param_name);
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = ['_' rlz '_' grid_label];
    % define filepaths
    nc_filepath_abs_sal = [fpath 'combined/' grid_type '/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param,...
        units,long_param_name);
end

%% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% start timing predictions
tStart = tic;

%% compute and save estimates for each month
parfor m = 1:length(TS.time)
    if strcmp(base_grid,'RG')
        % counter 
        cnt = m;
        % load dimensions
        TS = load_RG_dim([fpath '/Data/RG_CLIM/']);
        TS = replicate_dims(base_grid,TS,1);
        TS.longitude = convert_lon(TS.longitude);
        % get RG T and S
        TS.temperature = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Sal.nc'],'Salinity',...
            [1 1 1 m],[Inf Inf Inf 1]);
        % covert RG T and S to conservative temperature and absolute salinity
        pres_3d = repmat(permute(TS.Pressure,[3 2 1]),length(TS.Longitude),length(TS.Latitude),1);
        lon_3d = repmat(TS.Longitude,1,length(TS.Latitude),length(TS.Pressure));
        lat_3d = repmat(TS.Latitude',length(TS.Longitude),1,length(TS.Pressure));
        TS.salinity_abs = gsw_SA_from_SP(TS.salinity,pres_3d,convert_lon(lon_3d),lat_3d);
        TS.temperature_cns = gsw_CT_from_t(TS.salinity_abs,TS.temperature,pres_3d);
        % get time variables for just this timestep
        TS.Time = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time',m,1);
        date_temp = datevec(datenum(2004,1,1+double(TS.Time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply model
        apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
            base_grid,m,1,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_alg_dir,param,units,long_param_name);
    elseif strcmp(base_grid,'RFROM')
        cnt = m;
        % load dimensions
        TS = load_RFROM_dim([fpath '/Data/RFROM/']);
        TS = replicate_dims(base_grid,TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % counter 
            cnt = cnt + 1;
            % get RFROM T and S
            TS.temperature_cns = ncread([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity_abs = ncread([fpath '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % get time variables for just this timestep
            date_temp = datevec(datenum(TS.years(m),TS.months(m),15));
            date_temp0 = date_temp;
            date_temp0(:,2:3) = 1; % Jan. 1 of each year
            TS.year = date_temp(:,1);
            TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
            % transform day
            TS.day_sin = sin((2.*pi.*TS.day)/365.25);
            TS.day_cos = cos((2.*pi.*TS.day)/365.25);
            % apply model
            apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
                base_grid,m,w,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
                thresh,gobai_alg_dir,param,units,long_param_name);
        end
    else
        % counter 
        cnt = m;
        % define paths
        path2 = ['_Omon_' base_grid '_'];
        path3 = ['_'  rlz '_' grid_label];
        % define filepaths
        nc_filepath_abs_sal = [fpath 'combined/' grid_type '/abs_sal' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        nc_filepath_cns_tmp = [fpath 'combined/' grid_type '/cns_tmp' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        % load dimensions
        TS = load_model_dim(nc_filepath_abs_sal);
        TS = replicate_dims(base_grid,TS,1);
        % get practical salinity and in situ temperature from cmip model
        TS.salinity_abs = ncread(nc_filepath_abs_sal,'abs_sal',[1 1 1 m],[Inf Inf Inf 1]);
        TS.temperature_cns = ncread(nc_filepath_cns_tmp,'cns_tmp',[1 1 1 m],[Inf Inf Inf 1]);
        % get time variables for just this timestep
        TS.Time = ncread(nc_filepath_abs_sal,'time',m,1);
        date_temp = datevec(datenum(0,0,double(TS.Time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply model
        apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
            base_grid,m,1,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_alg_dir,param,units,long_param_name);
    end
end

% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf([alg_type ' Prediction (' num2str(start_year) ' to ' date_str(1:4) '): ']);

tElapsed = toc(tStart);
disp(['Elapsed time is ' num2str(tElapsed/60) ' minutes.'])

end

%% embedded functions

% for processing 3D grids and applying trained models to them
function apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
    base_grid,m,w,cnt,xdim,ydim,zdim,variables_TS,thresh,gobai_alg_dir,...
    param,units,long_param_name)
    
    % define file name
    filename = [gobai_alg_dir 'gobai-' param '.nc'];

    % convert to arrays
    TS_index = ~isnan(TS.temperature_cns);
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if ndims(TS.(vars{v})) == 3
            TS.([vars{v} '_array']) = TS.(vars{v})(TS_index);
            TS = rmfield(TS,vars{v});
        end
    end

    % replicate time variables as arrays
    vars = fieldnames(TS);
    for v = 1:length(vars)
        if length(TS.(vars{v})) == 1
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.temperature_cns_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute salinity, conservative temperature, potential density
    TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temperature_cns_array),num_clusters));
    probs_matrix = single(nan(length(TS.temperature_cns_array),num_clusters));

    % apply models for each cluster
    for c = 1:num_clusters

      % check for existence of model
      if isfile([alg_dir '/' alg_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        GMM_probs = ...
            load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/c' num2str(c) ...
            '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
        GMM_probs.probabilities_array = GMM_probs.GMM_cluster_probs(TS_index); % convert to array
        GMM_probs.probabilities_array(GMM_probs.probabilities_array < thresh) = NaN; % remove probabilities below thresh
        GMM_probs = rmfield(GMM_probs,{'GMM_cluster_probs'}); % remove 3D matrix
        probs_matrix(:,c) = GMM_probs.probabilities_array; % add to probability matrix

        % load model for this cluster
        alg_struct = load([alg_dir '/' alg_fnames{c}],alg_type);
    
        % predict data for each cluster
        if strcmp(alg_type,'FFNN')
            gobai_matrix(:,c) = ...
                run_FFNN(alg_struct.FFNN,TS,GMM_probs.probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS,thresh);
        elseif strcmp(alg_type,'RFR')
            gobai_matrix(:,c) = ...
                run_RFR(alg_struct.RFR,TS,GMM_probs.probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS,thresh);
        elseif strcmp(alg_type,'GBM')
            gobai_matrix(:,c) = ...
                run_GBM(alg_struct.GBM,TS,GMM_probs.probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS,thresh);
        end

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;

    %% Write monthly output
    ncwrite(filename,'time',datenum(TS.years(cnt),TS.months(cnt),15),cnt);
    ncwrite(filename,param,gobai_3d,[1 1 1 cnt]);

    % display information
    fprintf([alg_type ' Prediction (Month ' num2str(m) ', Week ' num2str(w) ')\n']);

end

% for creating netCDF file
function create_nc_file(TS,base_grid,xdim,ydim,zdim,gobai_alg_dir,...
    param,units,long_param_name)

% define file name
filename = [gobai_alg_dir 'gobai-' param '.nc'];

% create folder and file
if ~isfolder([pwd '/' gobai_alg_dir]); mkdir(gobai_alg_dir); end
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
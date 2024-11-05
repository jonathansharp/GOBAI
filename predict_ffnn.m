% predict_ffnn
%
% DESCRIPTION:
% This function uses trained FFNN models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/31/2024

function predict_ffnn(param,fpath,base_grid,file_date,float_file_ext,...
    num_clusters,variables,train_ratio,val_ratio,test_ratio,...
    thresh,numWorkers_predict,start_year,snap_date,varargin)

%% process date
date_str = num2str(snap_date);

%% directory base
dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
    float_file_ext;train_ratio;val_ratio;test_ratio});

%% process parameter name
[param1,param2,~,~,~,~,~,~,units,long_param_name] = param_name(param);

%% create directory and file names
ffnn_dir = [param1 '/Models/' dir_base];
ffnn_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    ffnn_fnames(c) = ...
        {['FFNN_' param2 '_C' num2str(c)]};
end
gobai_ffnn_dir = ...
    [param1 '/Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
    '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
    '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) '/'];

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% determine timesteps
% if strcmp(base_grid,'RG')
%     [~,timesteps] = load_RG_dim([fpath '/Data/RG_CLIM/']);
% elseif strcmp(base_grid,'RFROM')
%     [~,timesteps] = load_RFROM_dim([fpath '/Data/RFROM/']);
% end

%% predict property on grid

% set up parallel pool
% tic; parpool(numWorkers_predict); fprintf('Pool initiation:'); toc;

% start timing predictions
tic

% compute estimates for each month
m_last = (str2num(date_str(1:4))-2004)*12+str2num(date_str(5:6));
for m = 1:m_last
    if strcmp(base_grid,'RG')
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
        % apply ffnn model
        apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
            base_grid,m,1,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_ffnn_dir,param,units,long_param_name);
    elseif strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim([fpath '/Data/RFROM/']);
        TS = replicate_dims(base_grid,TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
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
            % apply ffnn model
            apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
                base_grid,m,w,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
                thresh,gobai_ffnn_dir,param,units,long_param_name);
        end
    else
        % define paths
        path2 = ['_Omon_' base_grid '_'];
        path3 = '_r1i1p1f1_gr';
        % define filepaths
        nc_filepath_abs_sal = [fpath 'combined/regridded/abs_sal' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        nc_filepath_cns_tmp = [fpath 'combined/regridded/cns_tmp' path2 ...
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
        % apply ffnn model
        apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
            base_grid,m,1,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_ffnn_dir,param,units,long_param_name);
    end
end

% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf(['FFNN Prediction (' num2str(start_year) ' to ' date_str(1:4) '): ']);
toc

end

%% embedded function for processing 3D grids and applying FFNN models
function apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
    base_grid,m,w,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir,...
    param,units,long_param_name)

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
      if isfile([ffnn_dir '/' ffnn_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        GMM_probs = ...
            load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/c' num2str(c) ...
            '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
        GMM_probs.probabilities_array = GMM_probs.GMM_cluster_probs(TS_index); % convert to array
        GMM_probs.probabilities_array(GMM_probs.probabilities_array < thresh) = NaN; % remove probabilities below thresh
        GMM_probs = rmfield(GMM_probs,{'GMM_cluster_probs'}); % remove 3D matrix
        probs_matrix(:,c) = GMM_probs.probabilities_array; % add to probability matrix

        % load model for this cluster
        alg = load([ffnn_dir '/' ffnn_fnames{c}],'FFNN');
    
        % predict data for each cluster
        gobai_matrix(:,c) = ...
            run_FFNN(alg.FFNN,TS,GMM_probs.probabilities_array,...
            true(size(TS.temperature_cns_array)),variables_TS,thresh);

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;
    
    % create folder and save monthly output
    if ~isfolder([pwd '/' gobai_ffnn_dir]); mkdir(gobai_ffnn_dir); end
    filename = [gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w) '.nc'];
    if isfile(filename); delete(filename); end % delete file if it exists
    % parameter
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        nccreate(filename,param,'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim},...
            'DataType','single','FillValue',NaN);
    else
        nccreate(filename,param,'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim},...
            'DataType','single','FillValue',NaN);
    end
    ncwrite(filename,param,gobai_3d);
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

    % display information
    fprintf(['FFNN Prediction (Month ' num2str(m) ', Week ' num2str(w) ')\n']);

end
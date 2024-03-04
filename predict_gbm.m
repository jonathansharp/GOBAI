% predict_gbm
%
% DESCRIPTION:
% This function uses trained GBM models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/31/2024

%% create directory and file names
gbm_dir = ['Models/' dir_base];
gbm_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    gbm_fnames(c) = ...
        {['GBM_oxygen_C' num2str(c)]};
end
gobai_gbm_dir = ...
    ['Data/GOBAI/' base_grid '/GBM/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numstumps) '/'];

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% determine timesteps
if strcmp(base_grid,'RG')
    [~,timesteps] = load_RG_dim([pwd '/Data/RG_CLIM/']);
elseif strcmp(base_grid,'RFROM')
    [~,timesteps] = load_RFROM_dim([pwd '/Data/RFROM/']);
end

%% predict property on grid

% start timing predictions
tic

% set up parallel pool
tic; parpool(12); fprintf('Pool initiation:'); toc;

% compute estimates for each month
m1 = (years_to_predict(1)-2004)*12+1;
m2 = (years_to_predict(end)-2004)*12+12;
parfor m = m1:m2
    if strcmp(base_grid,'RG')
        % load dimensions
        TS = load_RG_dim([pwd '/Data/RG_CLIM/']);
        TS = replicate_RG_dim(TS,1);
        TS.longitude = convert_lon(TS.longitude);
        % get RG T and S
        TS.temperature = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Sal.nc'],'Salinity',...
            [1 1 1 m],[Inf Inf Inf 1]);
        % get RG time variables
        TS.Time = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time',m,1);
        date_temp = datevec(datenum(2004,1,1+double(TS.Time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply gbm model
        apply_gbm_model(TS,num_clusters,gbm_dir,gbm_fnames,...
            base_grid,m,1,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_gbm_dir);
    elseif strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim([pwd '/Data/RFROM/']);
        TS = replicate_RFROM_dim(TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM T and S
            TS.temperature = ncread([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity = ncread([pwd '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            TS.time = ncread([pwd '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'time',w,1);
            % get RFROM time variables
            date_temp = datevec(datenum(1950,1,1+double(TS.time)));
            date_temp0 = date_temp;
            date_temp0(:,2:3) = 1; % Jan. 1 of each year
            TS.year = date_temp(:,1);
            TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
            % transform day
            TS.day_sin = sin((2.*pi.*TS.day)/365.25);
            TS.day_cos = cos((2.*pi.*TS.day)/365.25);
            % apply gbm model
            apply_gbm_model(TS,num_clusters,gbm_dir,gbm_fnames,...
            base_grid,m,w,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_gbm_dir);
        end
    end
end


% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf('GBM Prediction: ');
toc

% clean up


%% embedded function for processing 3D grids and applying GBM models
function apply_gbm_model(TS,num_clusters,gbm_dir,gbm_fnames,...
    base_grid,m,w,xdim,ydim,zdim,variables_TS,thresh,gobai_gbm_dir)

    % convert to arrays
    TS_index = ~isnan(TS.temperature);
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
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.temperature_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute salinity, conservative temperature, potential density
    TS.sigma_array = gsw_sigma0(TS.salinity_array,TS.temperature_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temperature_array),num_clusters));
    probs_matrix = single(nan(length(TS.temperature_array),num_clusters));

    % apply models for each cluster
    for c = 1:num_clusters

      % check for existence of model
      if isfile([gbm_dir '/' gbm_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        GMM_probs = ...
            load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/c' num2str(c) ...
            '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
        GMM_probs.probabilities_array = GMM_probs.GMM_cluster_probs(TS_index);
        GMM_probs.probabilities_array(GMM_probs.probabilities_array < thresh) = NaN; % remove probabilities below thresh
        GMM_probs = rmfield(GMM_probs,{'GMM_cluster_probs'});
        probs_matrix(:,c) = GMM_probs.probabilities_array;
    
        % load model for this cluster
        alg = load([gbm_dir '/' gbm_fnames{c}],'GBM');
    
        % predict data for each cluster
        gobai_matrix(:,c) = ...
            run_GBM(alg.GBM,TS,GMM_probs.probabilities_array,...
            true(size(TS.temperature_array)),variables_TS,thresh);
    
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
    if ~isfolder([pwd '/' gobai_gbm_dir]); mkdir(gobai_gbm_dir); end
    filename = [gobai_gbm_dir 'm' num2str(m) '_w' num2str(w) '.nc'];
    if isfile(filename); delete(filename); end % delete file if it exists
    % oxygen
    nccreate(filename,'o2','Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'o2',gobai_3d);
    ncwriteatt(filename,'o2','units','umol/kg');
    ncwriteatt(filename,'o2','long_name','Dissolved Oxygen Amount Content');
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
    % pressure
    nccreate(filename,'pres','Dimensions',{'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'pres',TS.Pressure);
    ncwriteatt(filename,'pres','units','decibars');
    ncwriteatt(filename,'pres','axis','Z');
    ncwriteatt(filename,'pres','long_name','pressure');
    ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');

%     parsave([gobai_gbm_dir 'm' num2str(m) '_w' num2str(w)],...
%         gobai_3d,'gobai',w,'w',m,'m');

end

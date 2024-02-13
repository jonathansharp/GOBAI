% predict_ffnn (hackathon test)
%
% DESCRIPTION:
% This function uses trained FFNN models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/12/2024

%% load configuration parameters
gobai_o2_initiate;
load_standard_config_files;
load('Config/base_config_RFROM.mat');
load('Config/predict_years_config_04.mat');
dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
    float_file_ext;train_ratio;val_ratio;test_ratio});
fpath = '/raid';

%% create directory and file names
ffnn_dir = ['Models/' dir_base];
ffnn_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    ffnn_fnames(c) = ...
        {['FFNN_oxygen_C' num2str(c)]};
end
gobai_ffnn_dir = ...
    ['Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% determine timesteps
if strcmp(base_grid,'RG')
    [~,timesteps] = load_RG_dim([fpath '/Data/RG_CLIM/']);
elseif strcmp(base_grid,'RFROM')
    [~,timesteps] = load_RFROM_dim([fpath '/Data/RFROM/']);
end

%% predict property on grid

% start timing predictions
tic

% set up parallel pool
% tic; parpool(12); fprintf('Pool initiation:'); toc;

% compute estimates for each month
m1 = (years_to_predict(1)-2004)*12+1;
m2 = (years_to_predict(end)-2004)*12+12;
% parfor m = m1:m2
for m = m1:m2
    if strcmp(base_grid,'RG')
        % load dimensions
        TS = load_RG_dim([fpath '/Data/RG_CLIM/']);
        TS = replicate_RG_dim(TS,1);
        TS.longitude = convert_lon(TS.longitude);
        % get RG T and S
        TS.temperature = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([fpath '/Data/RG_CLIM/RG_Climatology_Sal.nc'],'Salinity',...
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
        % apply ffnn model
        apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
            base_grid,m,1,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_ffnn_dir);
    elseif strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim([fpath '/Data/RFROM/']);
        TS = replicate_RFROM_dim(TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM T and S
            TS.temperature = ncread([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity = ncread([fpath '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % get RFROM time variables
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
            thresh,gobai_ffnn_dir);
        end
    end
end


% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf(['FFNN Prediction (' num2str(years_to_predict(1)) ' to ' num2str(years_to_predict(end)) '): ']);
toc

% clean up

%% determine seasonal averages

%% determine annual averages

%% determine long-term average

%% embedded function for processing 3D grids and applying FFNN models
function apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
    base_grid,m,w,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir)

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
    TS.salinity_abs_array = gsw_SA_from_SP(TS.salinity_array,TS.pressure_array,...
        TS.longitude_array,TS.latitude_array);
    TS.temperature_cns_array = gsw_CT_from_t(TS.salinity_abs_array,...
        TS.temperature_array,TS.pressure_array);
    TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temperature_array),num_clusters));
    probs_matrix = single(nan(length(TS.temperature_array),num_clusters));

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
    if ~isfolder([fpath '/' gobai_ffnn_dir]); mkdir(gobai_ffnn_dir); end
    filename = [gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w) '.nc'];
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
    
%     parsave([gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w)],...
%         gobai_3d,'gobai',w,'w',m,'m');

end
% train_ffnn
%
% DESCRIPTION:
% This function uses the combined dataset to train
% neural networks in GMM clusters and predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/3/2024

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' base_grid '_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% remove float data for GLODAP only test
if glodap_only
    glodap_idx = all_data.platform < 10^6;
    vars = fieldnames(all_data);
    for v = 1:length(vars)
        all_data.(vars{v}) = all_data.(vars{v})(glodap_idx);
    end
end
clear glodap_idx vars v

%% create directory and file names
ffnn_dir = ['Models/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
ffnn_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    ffnn_fnames(c) = ...
        {['FFNN_oxygen_C' num2str(c)]};
end
gobai_ffnn_dir = ...
    ['Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];

%% fit FFNNs using all data

% start timing training
tic

% define model parameters
nodes1 = [5 10 15];
nodes2 = [15 10 5];

% fit models for each cluster
for c = 1:num_clusters

  % check for data in cluster
  if any(all_data_clusters.clusters == c)
    
    % start timing fit
    tic

    % fit model for each cluster
    FFNN = ...
        fit_FFNN('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
        true(size(all_data.platform)),variables,nodes1,nodes2,...
        train_ratio,val_ratio,test_ratio,thresh);

    % save model for each cluster
    if ~isfolder([pwd '/' ffnn_dir]); mkdir(ffnn_dir);end
    save([ffnn_dir '/' ffnn_fnames{c}],'FFNN','-v7.3');

    % clean up
    clear FFNN

    % stop timing fit
    fprintf(['Train FFNN - Cluster #' num2str(c) ': ']);
    toc

  else

    % stop timing fit
    fprintf(['Train FFNN - Cluster #' num2str(c) ': N/A']);
    [~]=toc;

  end

end

% clean up
clear all_data all_data_clusters

% end parallel session
delete(gcp('nocreate'));

% stop timing training
fprintf('FFNN Training: ');
toc

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% define file path and load dimensional variables
if strcmp(base_grid,'RG')

    % file path
    fpath = [pwd '/Data/RG_CLIM/'];

    % load dimensions
    lon = ncread([fpath 'RG_Climatology_Temp.nc'],'Longitude');
    lat = ncread([fpath 'RG_Climatology_Temp.nc'],'Latitude');
    pres = ncread([fpath 'RG_Climatology_Temp.nc'],'Pressure');
    lon = convert_lon(lon);

    % compute dimensions
    xdim = length(lon);
    ydim = length(lat);
    zdim = length(pres);    

    % get timesteps
    time_temp = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time');
    timesteps = length(time_temp);
    clear time_temp

elseif strcmp(base_grid,'RFROM')

    % file path
    fpath = [pwd '/Data/RFROM/'];

    % load RFROM climatological temp and salinity
    lon = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'longitude');
    lat = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'latitude');
    pres = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'mean_pressure');
    
    % compute dimensions
    xdim = length(lon);
    ydim = length(lat);
    zdim = length(pres);

    % get timesteps
    files = dir([fpath 'RFROM_TEMP_v0.1/*.nc']);
    timesteps = length(files);
    clear files

end

% expand latitude, longitude, and depth
longitude = repmat(lon,1,ydim,zdim);
lon_cos_1 = repmat(cosd(lon-20),1,ydim,zdim);
lon_cos_2 = repmat(cosd(lon-110),1,ydim,zdim);
latitude = repmat(lat',xdim,1,zdim);
pressure = repmat(permute(pres,[3 2 1]),xdim,ydim,1);

%% process time for RFROM
if strcmp(base_grid,'RFROM')
    files = dir([fpath 'RFROM_TEMP_v0.1/*.nc']);
    years = nan(length(files),1);
    months = nan(length(files),1);
    for n = 1:length(files)
        date = cell2mat(extractBetween(files(n).name,'STABLE_','.nc'));
        years(n) = str2double(date(1:4));
        months(n) = str2double(date(6:7));
    end
    idx = years < 2004;
    years(idx) = [];
    months(idx) = [];
end

%% predict property on grid

% start timing predictions
tic

% set up parallel pool
%parpool;

% compute estimates for each month
for m = 1:timesteps
    % populate structure
    TS = struct();
    TS.longitude = longitude;
    TS.lon_cos_1 = lon_cos_1;
    TS.lon_cos_2 = lon_cos_2;
    TS.latitude = latitude;
    TS.pressure = pressure;
    % load T/S grid
    if strcmp(base_grid,'RG')
        % get RG T and S
        TS.temp = ncread([fpath 'RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
        TS.sal = ncread([fpath 'RG_Climatology_Sal.nc'],'Salinity',...
            [1 1 1 m],[Inf Inf Inf 1]);
        % get RG time variables
        TS.time = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time',...
            m,1);
        date_temp = datevec(datenum(2004,1,1+double(TS.time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply ffnn model
        apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
            base_grid,m,1,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir);
    elseif strcmp(base_grid,'RFROM')
        % determine number of weeks in file
        nc_atts = ncinfo([fpath 'RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(years(m)) '_' sprintf('%02d',months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM T and S
            TS.temp = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.sal = ncread([fpath 'RFROM_SAL_STABLE_CLIM.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % get RFROM time variables
            date_temp = datevec(datenum(years(m),months(m),15));
            date_temp0 = date_temp;
            date_temp0(:,2:3) = 1; % Jan. 1 of each year
            TS.year = date_temp(:,1);
            TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
            % transform day
            TS.day_sin = sin((2.*pi.*TS.day)/365.25);
            TS.day_cos = cos((2.*pi.*TS.day)/365.25);
            % apply ffnn model
            apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
            base_grid,m,1,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir);
        end
    end
end


% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf('FFNN Prediction: ');
toc

% clean up


%% embedded function for processing 3D grids and applying FFNN models
function apply_ffnn_model(TS,num_clusters,ffnn_dir,ffnn_fnames,...
    base_grid,m,w,xdim,ydim,zdim,variables_TS,thresh,gobai_ffnn_dir)

    % convert to arrays
    TS_index = ~isnan(TS.temp);
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
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.temp_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute salinity, conservative temperature, potential density
    TS.salinity_abs_array = gsw_SA_from_SP(TS.sal_array,TS.pressure_array,...
        TS.longitude_array,TS.latitude_array);
    TS.temperature_cns_array = gsw_CT_from_t(TS.salinity_abs_array,...
        TS.temp_array,TS.pressure_array);
    TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temp_array),num_clusters));
    probs_matrix = single(nan(length(TS.temp_array),num_clusters));

    % apply models for each cluster
    for c = 1:2%num_clusters

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
            true(size(TS.temp_array)),variables_TS,thresh);

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;
    
    % save monthly output
    if ~isfolder([pwd '/' gobai_ffnn_dir]); mkdir(gobai_ffnn_dir); end
    parsave([gobai_ffnn_dir 'm' num2str(m) '_w' num2str(w)],...
        gobai_3d,'gobai',w,'w',m,'m');

end

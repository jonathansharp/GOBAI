% predict_rfr
%
% DESCRIPTION:
% This function uses trained RFR models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/28/2024

%% initiate profile
profile on

%% create directory and file names
rfr_dir = ['Models/' dir_base];
rfr_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    rfr_fnames(c) = ...
        {['RFR_oxygen_C' num2str(c)]};
end
gobai_rfr_dir = ...
    ['Data/GOBAI/' base_grid '/RFR/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];
if ~isfolder([pwd '/' gobai_rfr_dir]); mkdir(gobai_rfr_dir); end

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% predict property on grid

% start timing predictions
tic

% % set up parallel pool
tic; parpool(12); fprintf('Pool initiation:'); toc;
[X,message_ID] = lastwarn()

% parallelize across clusters
parfor c = 1:num_clusters

    % check for existence of model
    if isfile([rfr_dir '/' rfr_fnames{c} '.mat'])

        % load model for this cluster
        alg = load([rfr_dir '/' rfr_fnames{c}],'RFR');
    
        % compute estimates for each month
        m1 = (years_to_predict(1)-2004)*12+1;
        m2 = (years_to_predict(end)-2004)*12+12;
        for m = m1:m2
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
                % apply rfr model
                apply_rfr_model(TS,num_clusters,...
                    base_grid,m,1,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
                    thresh,gobai_rfr_dir,c);
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
                    % get RFROM time variables
                    date_temp = datevec(datenum(TS.years(m),TS.months(m),15));
                    date_temp0 = date_temp;
                    date_temp0(:,2:3) = 1; % Jan. 1 of each year
                    TS.year = date_temp(:,1);
                    TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
                    % transform day
                    TS.day_sin = sin((2.*pi.*TS.day)/365.25);
                    TS.day_cos = cos((2.*pi.*TS.day)/365.25);
                    % apply rfr model
                    apply_rfr_model(alg,TS,num_clusters,...
                        base_grid,m,w,variables_TS,thresh,gobai_rfr_dir,c);
                    % print status
                    fprintf(['Month #' num2str(m) ', Week #' num2str(w) ' Done']);
                end
            end
        end

    else



    end

end

% end parallel session
delete(gcp('nocreate'));

% stop timing predictions
fprintf('RFR Prediction: ');
toc

% clean up

%% reform 3D grids and save as NetCDFs

% start timing predictions
tic

% set up parallel pool
tic; parpool(12); fprintf('Pool initiation:'); toc;

% parallelize across months
m1 = (years_to_predict(1)-2004)*12+1;
m2 = (years_to_predict(end)-2004)*12+12;
for m = m1:m2

    if strcmp(base_grid,'RG')
        % load dimensions
        TS = load_RG_dim([pwd '/Data/RG_CLIM/']);
        TS = replicate_RG_dim(TS,1);
        TS.longitude = convert_lon(TS.longitude);
        % get RG T and S
        TS.temperature = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
    elseif strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim([pwd '/Data/RFROM/']);
        TS = replicate_RFROM_dim(TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM
            TS.temperature = ncread([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
        end
    end

    % pre-allocate
    TS_index = ~isnan(TS.temperature);
    gobai_matrix = nan(sum(TS_index(:)),1);
    probs_matrix = nan(sum(TS_index(:)),1);

    % for each cluster
    for c = 1:num_clusters
        % check for existence of model
        if isfile([rfr_dir '/' rfr_fnames{c} '.mat'])
            fname = [gobai_rfr_dir 'probs_c' num2str(c) '_m' num2str(m) '_w' num2str(w) '.nc'];
            probs_matrix(:,c) = ncread(fname,'probs');
            fname = [gobai_rfr_dir 'preds_c' num2str(c) '_m' num2str(m) '_w' num2str(w) '.nc'];
            gobai_matrix(:,c) = ncread(fname,'preds');
        end
    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));

    % convert back to 3D grid
    gobai_3d = nan(TS.xdim,TS.ydim,TS.zdim);
    gobai_3d(TS_index) = gobai_array;

    % create folder and save as NetCDF
    if ~isfolder([pwd '/' gobai_rfr_dir]); mkdir(gobai_rfr_dir); end
    filename = [gobai_rfr_dir 'm' num2str(m) '_w' num2str(w) '.nc'];
    if isfile(filename); delete(filename); end % delete file if it exists
    % oxygen
    nccreate(filename,'o2','Dimensions',{'lon',TS.xdim,'lat',TS.ydim,'pres',TS.zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'o2',gobai_3d);
    ncwriteatt(filename,'o2','units','umol/kg');
    ncwriteatt(filename,'o2','long_name','Dissolved Oxygen Amount Content');
    % longitude
    nccreate(filename,'lon','Dimensions',{'lon',TS.xdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lon',TS.Longitude);
    ncwriteatt(filename,'lon','units','degrees_east');
    ncwriteatt(filename,'lon','axis','X');
    ncwriteatt(filename,'lon','long_name','longitude');
    ncwriteatt(filename,'lon','_CoordinateAxisType','Lon');
    % latitude
    nccreate(filename,'lat','Dimensions',{'lat',TS.ydim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'lat',TS.Latitude);
    ncwriteatt(filename,'lat','units','degrees_north');
    ncwriteatt(filename,'lat','axis','Y');
    ncwriteatt(filename,'lat','long_name','latitude');
    ncwriteatt(filename,'lat','_CoordinateAxisType','Lat');
    % pressure
    nccreate(filename,'pres','Dimensions',{'pres',TS.zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'pres',TS.Pressure);
    ncwriteatt(filename,'pres','units','decibars');
    ncwriteatt(filename,'pres','axis','Z');
    ncwriteatt(filename,'pres','long_name','pressure');
    ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');

    % print status
    fprintf(['Month #' num2str(m) ', Week #' num2str(w) ' Done']);

end

% end parallel session
delete(gcp('nocreate'));

% stop timing re-formation
fprintf('Grid Formaton: ');
toc

%% end and save profile
p=profile('info');
profsave(p,'profiles/train_rfr_new')
profile off

%% embedded function for processing 3D grids and applying RFR models
function apply_rfr_model(alg,TS,num_clusters,...
    base_grid,m,w,variables_TS,thresh,gobai_rfr_dir,c)

    % convert grids to arrays
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

    % load GMM cluster probabilities for this cluster and month, and convert to array
    GMM_probs = load(['Data/GMM_' base_grid '_' num2str(num_clusters) ...
        '/c' num2str(c) '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs');
    probabilities_array = GMM_probs.GMM_cluster_probs(TS_index);
    % remove probabilities below thresh
    probabilities_array(probabilities_array < thresh) = NaN;
    clear GMM_probs
    
    % predict data for each cluster
    gobai_array = ...
        run_RFR(alg.RFR,TS,probabilities_array,...
        true(size(TS.temperature_array)),variables_TS,thresh);
    
    % save probabilities
    fname = [gobai_rfr_dir 'probs_c' num2str(c) '_m' num2str(m) '_w' num2str(w) '.nc'];
    if exist(fname,'file'); delete(fname); end
    nccreate(fname,'probs','Dimensions',{'r' length(probabilities_array) 'c' 1});
    ncwrite(fname,'probs',probabilities_array);
    
    % save predictions
    fname = [gobai_rfr_dir 'preds_c' num2str(c) '_m' num2str(m) '_w' num2str(w) '.nc'];
    if exist(fname,'file'); delete(fname); end
    nccreate(fname,'preds','Dimensions',{'r' length(gobai_array) 'c' 1});
    ncwrite(fname,'preds',gobai_array);
    
end

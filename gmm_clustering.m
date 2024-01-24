% gmm_clustering
%
% DESCRIPTION:
% This function uses gridded temperature and salinity with Gaussian Mixture
% Modelling to formulate global clusters of similar environmental
% variability within which to train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 12/1/2023

%% load temperature and salinity data
if strcmp(base_grid,'RG')
    % file path
    fpath = [pwd '/Data/RG_CLIM/'];
    % load temperature and salinity
    TS.Longitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Longitude');
    TS.Latitude = ncread([fpath 'RG_Climatology_Temp.nc'],'Latitude');
    TS.Pressure = ncread([fpath 'RG_Climatology_Temp.nc'],'Pressure');
    TS.Temperature = ncread([fpath 'RG_Climatology_Temp.nc'],'Temperature');
    TS.Salinity = ncread([fpath 'RG_Climatology_Sal.nc'],'Salinity');
    % compute dimensions
    xdim = length(TS.Longitude);
    ydim = length(TS.Latitude);
    zdim = length(TS.Pressure);
    % calculate climatological mean temperature and salinity
    TS.temp_clim = single(nan(xdim,ydim,zdim,12));
    TS.sal_clim = single(nan(xdim,ydim,zdim,12));
    for m = 1:12
        TS.temp_clim(:,:,:,m) = mean(TS.Temperature(:,:,:,m:12:end),4,'omitnan');
        TS.sal_clim(:,:,:,m) = mean(TS.Salinity(:,:,:,m:12:end),4,'omitnan');
    end
    % clean up
    TS = rmfield(TS,{'Temperature' 'Salinity'});
    % determine number of monthly timesteps
    TS.Time = ncread('Data/RG_CLIM/RG_Climatology_Temp.nc','Time');
    timesteps = length(TS.Time);
elseif strcmp(base_grid,'RFROM')
    % file path
    fpath = [pwd '/Data/RFROM/'];
    % load RFROM climatological temp and salinity
    TS.Longitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'longitude');
    TS.Latitude = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'latitude');
    TS.Pressure = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'mean_pressure');
    TS.temp_clim = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],'ocean_temperature');
    TS.sal_clim = ncread([fpath 'RFROM_SAL_STABLE_CLIM.nc'],'ocean_salinity');
    % compute dimensions
    xdim = length(TS.Longitude);
    ydim = length(TS.Latitude);
    zdim = length(TS.Pressure);
    % determine number of monthly timesteps
    timesteps = length(dir([fpath 'RFROM_TEMP_v0.1/*.nc']));
end

%% compute long-term means of T and S (for RFROM cause it's big)
if strcmp(base_grid,'RFROM')
    TS.temp_clim = mean(TS.temp_clim,4,'omitnan');
    TS.sal_clim = mean(TS.sal_clim,4,'omitnan');
end

%% expand latitude, longitude, and depth
% just annual mean for RFROM but monthly means for RG
if strcmp(base_grid,'RG')
    tdim = 12;
elseif strcmp(base_grid,'RFROM')
    tdim = 1;
end
TS.lon_cos_1_3D = repmat(cosd(TS.Longitude-20),1,ydim,zdim,tdim);
TS.lon_cos_2_3D = repmat(cosd(TS.Longitude-110),1,ydim,zdim,tdim);
TS.latitude_3D = repmat(TS.Latitude',xdim,1,zdim,tdim);
TS.pressure_3D = repmat(permute(TS.Pressure,[3 2 1]),xdim,ydim,1,tdim);

%% fit GMM from climatological mean temperature and salinity
% tic
% % transform to normalized arrays
% idx = ~isnan(TS.temp_clim) & ~isnan(TS.sal_clim);
% predictor_matrix = [];
% for v = 1:length(clust_vars_grid_1)
%     predictor_matrix = [predictor_matrix TS.(clust_vars_grid_1{v})(idx)];
% end
% [X_norm,C,S] = normalize(predictor_matrix);
% TS = rmfield(TS,{'temp_clim' 'sal_clim' 'lon_cos_1_3D' 'lon_cos_2_3D' 'latitude_3D' 'pressure_3D'});
% % fit GMM
% gmm = fitgmdist(X_norm,num_clusters,...
%     'CovarianceType','full',...
%     'SharedCovariance',true,'Replicates',10);
% % save GMM model
% if ~isfolder(['Data/GMM_' base_grid '_' num2str(num_clusters)])
%     mkdir(['Data/GMM_' base_grid '_' num2str(num_clusters)]);
% end
% save(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model'],...
%     'gmm','num_clusters','C','S','-v7.3');
% clear gmm C S
% toc

%% process time
if strcmp(base_grid,'RG')
    dates = datevec(datenum(2004,1,1+double(TS.Time)));
    years = dates(:,1);
    months = dates(:,2);
    clear dates
elseif strcmp(base_grid,'RFROM')
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

%% assign grid cells and probabilities to clusters

% start timing cluster assignment
tic

% load GMM model
load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model'],'gmm','C','S');

% remove dimensional variables from structure
lon_cos_1_3D = TS.lon_cos_1_3D;
lon_cos_2_3D = TS.lon_cos_2_3D;
latitude_3D = TS.latitude_3D;
pressure_3D = TS.pressure_3D;
clear TS

% set up parallel pool
parpool;
% for each month
parfor m = 1:timesteps
    % load T/S grid
    if strcmp(base_grid,'RG')
        % get RG T/S
        temp_3D = ncread([fpath 'RG_Climatology_Temp.nc'],'Temperature',...
            [1 1 1 m],[Inf Inf Inf 1]);
        sal_3D = ncread([fpath 'RG_Climatology_Sal.nc'],'Salinity',...
            [1 1 1 m],[Inf Inf Inf 1]);
        % transform to normalized arrays
        idx = ~isnan(temp_3D) & ~isnan(sal_3D);
        predictor_matrix = [];
        for v = 1:length(clust_vars_grid_2)
            eval(['predictor_matrix = [predictor_matrix ' clust_vars_grid_2{v}']');
        end
        [X_norm,C,S] = normalize(predictor_matrix);
        X_norm = normalize([temp(idx) sal(idx) pres(idx) ...
            cos_lon_1(idx) cos_lon_2(idx) lat(idx)],...
            'Center',C,'Scale',S);
        % assign data points to clusters
        assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,xdim,...
            ydim,zdim,m,1,base_grid);
    else
        % determine number of weeks in file
        nc_atts = ncinfo([fpath 'RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(years(m)) '_' sprintf('%02d',months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM T and S
            temp_3D = ncread([fpath 'RFROM_TEMP_STABLE_CLIM.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            sal_3D = ncread([fpath 'RFROM_SAL_STABLE_CLIM.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % transform to normalized arrays
            idx = ~isnan(temp) & ~isnan(sal);
            X_norm = normalize([temp(idx) sal(idx) pres(idx)],...
                'Center',C,'Scale',S);
            % assign data points to clusters
            assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,...
                xdim,ydim,zdim,m,w,base_grid);
        end
    end
end
% end parallel session
delete(gcp('nocreate'));

toc

%% concatenate clusters and probabilities
% % load and concatenate clusters
% GMM.clusters = uint8([]);
% GMM.longitude = TS.Longitude;
% GMM.latitude = TS.Latitude;
% GMM.pressure = TS.Pressure;
% GMM.time = TS.Time;
% for m = 1:length(TS.Time)
%     load(['Data/GMM_' num2str(num_clusters) '/m' num2str(m)],'GMM_monthly');
%     GMM.clusters = cat(4,GMM.clusters,uint8(GMM_monthly));
%     clear GMM_monthly
% end
% if ~isfolder('Data'); mkdir('Data'); end
% save(['Data/GMM_' num2str(num_clusters) '/GMM'],'GMM','-v7.3');
% clear GMM
% % load and concatenate cluster probabilities
% for c = 1:num_clusters
%     GMM_probs.probabilities = single([]);
%     GMM_probs.longitude = TS.Longitude;
%     GMM_probs.latitude = TS.Latitude;
%     GMM_probs.pressure = TS.Pressure;
%     GMM_probs.time = TS.Time;
%     for m = 1:length(TS.Time)
%         load(['Data/GMM_' num2str(num_clusters) '/c' num2str(c) ...
%             '/m' num2str(m)],'GMM_monthly_probs');
%         GMM_probs.probabilities = ...
%             cat(4,GMM_probs.probabilities,single(GMM_monthly_probs));
%         clear GMM_monthly_probs
%     end
%     if ~isfolder(['Data/GMM_' num2str(num_clusters)])
%         mkdir(['Data/GMM_' num2str(num_clusters)]);
%     end
%     save(['Data/GMM_' num2str(num_clusters) '/GMM_c' num2str(c)],'GMM_probs','-v7.3');
%     clear GMM_probs
% end


% old time info
% 1.8 hours for ten clusters and ten replicates on chinook? (9/8/23)
% 2.2 hours for twenty clusters and twenty replicates on Hercules (9/8/23)

% embedded function
function assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,xdim,ydim,zdim,m,w,base_grid)
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters
    GMM_clusters = nan(xdim,ydim,zdim);
    GMM_clusters(idx) = clusters;
    % save clusters
    folder_name = [pwd '/Data/GMM_' base_grid '_' num2str(num_clusters)];
    if ~isfolder(folder_name); mkdir(folder_name); end
    save([folder_name '/m' num2str(m) '_w' num2str(w)],'GMM_clusters','-v7.3');
    clear GMM_clusters
    % fill 3D probabilities
    for c = 1:num_clusters
        GMM_cluster_probs = nan(xdim,ydim,zdim);
        GMM_cluster_probs(idx) = p(:,c);
        % save probabilities
        folder_name = [pwd '/Data/GMM_' base_grid '_' num2str(num_clusters) '/c' num2str(c)];
        if ~isfolder(folder_name); mkdir(folder_name); end
        save([folder_name '/m' num2str(m) '_w' num2str(w)],'GMM_cluster_probs','-v7.3');
        clear GMM_cluster_probs
    end
end

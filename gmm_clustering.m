% gmm_clustering_rg
%
% DESCRIPTION:
% This function uses gridded temperature and salinity with Gaussian Mixture
% Modelling to formulate global clusters of similar environmental
% variability within which to train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/17/2023

tic

%% load temperature and salinity data
if strcmp(base_grid,'RG')
    % load temperature and salinity
    TS = netcdfreader('Data/RG_CLIM/RG_Climatology_Temp.nc');
    TS.Salinity = ncread('Data/RG_CLIM/RG_Climatology_Sal.nc','Salinity');
    % compute dimensions
    xdim = length(TS.Longitude);
    ydim = length(TS.Latitude);
    zdim = length(TS.Pressure);
    tdim = length(TS.Time);
    % calculate climatological mean temperature and salinity
    TS.temp_clim = single(nan(xdim,ydim,zdim,12));
    TS.sal_clim = single(nan(xdim,ydim,zdim,12));
    for m = 1:12
        TS.temp_clim(:,:,:,m) = mean(TS.Temperature(:,:,:,m:12:end),4,'omitnan');
        TS.sal_clim(:,:,:,m) = mean(TS.Salinity(:,:,:,m:12:end),4,'omitnan');
    end
    % clean up
    TS = rmfield(TS,{'Temperature' 'Salinity'});
elseif strcmp(base_grid,'RFROM')
    % load RFROM climatological temp and salinity
    
end

%% expand latitude, longitude, and depth
TS.lon_cos_3D = repmat(cosd(TS.Longitude-20),1,ydim,zdim,12);
TS.latitude_3D = repmat(TS.Latitude',xdim,1,zdim,12);
TS.pressure_3D = repmat(permute(TS.Pressure,[3 2 1]),xdim,ydim,1,12);

%% fit GMM from climatological mean temperature and salinity
% transform to normalized arrays
idx = ~isnan(TS.temp_clim) & ~isnan(TS.sal_clim);
[X_norm,C,S] = normalize([TS.temp_clim(idx) TS.sal_clim(idx) TS.lon_cos_3D(idx)...
    TS.latitude_3D(idx) TS.pressure_3D(idx)]);
TS = rmfield(TS,{'temp_clim' 'sal_clim' 'lon_cos_3D' 'latitude_3D' 'pressure_3D'});
% fit GMM
gmm = fitgmdist(X_norm,num_clusters,...
    'CovarianceType','full',...
    'SharedCovariance',true,'Replicates',20);
% save GMM model
if ~isfolder(['Data/GMM_' num2str(num_clusters)])
    mkdir(['Data/GMM_' num2str(num_clusters)]);
end
save(['Data/GMM_' num2str(num_clusters) '/model'],...
    'gmm','num_clusters','C','S','-v7.3');
clear gmm C S

%% process time
start_year = 2004;
end_year = floor(snap_date/1e2);
if strcmp(base_grid,'RG')
    end_month = mod(snap_date,1e2);
    years = repelem(start_year:end_year,12)';
    months = repmat(1:12,1,length(years)/12)';
    years = years(1:end-(12-end_month));
    months = months(1:end-(12-end_month));
elseif strcmp(base_grid,'RFROM')
    end_month = mod(snap_date,1e2);
    years = repelem(start_year:end_year,12)';
    months = repmat(1:12,1,length(years)/12)';
    years = years(1:end-(12-end_month));
    months = months(1:end-(12-end_month));
end
clear start_year end_year end_month

%% assign grid cells and probabilities to clusters
% load GMM model
load(['Data/GMM_' num2str(num_clusters) '/model'],'gmm','C','S');
% for each timestep
cos_lon = repmat(cosd(TS.Longitude-20),1,length(TS.Latitude),length(TS.Pressure));
lat = repmat(TS.Latitude',length(TS.Longitude),1,length(TS.Pressure));
pres = repmat(permute(TS.Pressure,[3 2 1]),length(TS.Longitude),length(TS.Latitude),1);
for m = 1:length(TS.Time)
    % transform to normalized arrays
    temp = TS.Temperature(:,:,:,m);
    sal = TS.Salinity(:,:,:,m);
    idx = ~isnan(temp) & ~isnan(sal);
    X_norm = normalize([temp(idx) sal(idx) cos_lon(idx) lat(idx) pres(idx)],...
        'Center',C,'Scale',S);
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters
    GMM_monthly = nan(xdim,ydim,zdim);
    GMM_monthly(idx) = clusters;
    % save clusters
    if ~isfolder(['Data/GMM_' num2str(num_clusters)])
        mkdir(['Data/GMM_' num2str(num_clusters)]);
    end
    save(['Data/GMM_' num2str(num_clusters) '/m' num2str(m)],'GMM_monthly','-v7.3');
    clear GMM_monthly
    % fill 3D probabilities
    for c = 1:num_clusters
        GMM_monthly_probs = nan(xdim,ydim,zdim);
        GMM_monthly_probs(idx) = p(:,c);
        % save probabilities
        if ~isfolder(['Data/GMM_' num2str(num_clusters) '/c' num2str(c)])
            mkdir(['Data/GMM_' num2str(num_clusters) '/c' num2str(c)]);
        end
        save(['Data/GMM_' num2str(num_clusters) '/c' num2str(c) '/m' ...
            num2str(m)],'GMM_monthly_probs','-v7.3');
        clear GMM_monthly_probs
    end
end
% clean up
clear cos_lon lat pres

%% concatenate clusters and probabilities
% load and concatenate clusters
GMM.clusters = uint8([]);
GMM.longitude = TS.Longitude;
GMM.latitude = TS.Latitude;
GMM.pressure = TS.Pressure;
GMM.time = TS.Time;
for m = 1:length(TS.Time)
    load(['Data/GMM_' num2str(num_clusters) '/m' num2str(m)],'GMM_monthly');
    GMM.clusters = cat(4,GMM.clusters,uint8(GMM_monthly));
    clear GMM_monthly
end
if ~isfolder('Data'); mkdir('Data'); end
save(['Data/GMM_' num2str(num_clusters) '/GMM'],'GMM','-v7.3');
clear GMM
% load and concatenate cluster probabilities
for c = 1:num_clusters
    GMM_probs.probabilities = single([]);
    GMM_probs.longitude = TS.Longitude;
    GMM_probs.latitude = TS.Latitude;
    GMM_probs.pressure = TS.Pressure;
    GMM_probs.time = TS.Time;
    for m = 1:length(TS.Time)
        load(['Data/GMM_' num2str(num_clusters) '/c' num2str(c) ...
            '/m' num2str(m)],'GMM_monthly_probs');
        GMM_probs.probabilities = ...
            cat(4,GMM_probs.probabilities,single(GMM_monthly_probs));
        clear GMM_monthly_probs
    end
    if ~isfolder(['Data/GMM_' num2str(num_clusters)])
        mkdir(['Data/GMM_' num2str(num_clusters)]);
    end
    save(['Data/GMM_' num2str(num_clusters) '/GMM_c' num2str(c)],'GMM_probs','-v7.3');
    clear GMM_probs
end

% clean up
clear TS

toc
% 1.8 hours for ten clusters and ten replicates on chinook? (9/8/23)
% 2.2 hours for twenty clusters and twenty replicates on Hercules (9/8/23)

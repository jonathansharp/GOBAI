% gmm_clustering_rfrom
%
% DESCRIPTION:
% This function uses gridded temperature and salinity from RFROM maps
% (Lyman and Johnson, in prep) with Gaussian Mixture Modelling to formulate
% global clusters of similar environmental variability within which to
% train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 10/31/2023

tic

%% load temperature and salinity data
% import_RFROM
path = '/raid/Data/RFROM/';
load('/raid/Data/RFROM/RG_200401_202307.mat');

%% compute dimensions
xdim = length(RG.longitude);
ydim = length(RG.latitude);
zdim = length(RG.pressure);
tdim = length(RG.time);

%% calculate climatological mean temperature and salinity
RG.temp_clim = single(nan(xdim,ydim,zdim,12));
RG.sal_clim = single(nan(xdim,ydim,zdim,12));
for m = 1:12
    RG.temp_clim(:,:,:,m) = mean(RG.temp(:,:,:,m:12:end),4,'omitnan');
    RG.sal_clim(:,:,:,m) = mean(RG.sal(:,:,:,m:12:end),4,'omitnan');
end

%% expand latitude, longitude, and depth
RG.lon_cos_3D = repmat(cosd(RG.longitude-20),1,ydim,zdim,12);
RG.latitude_3D = repmat(RG.latitude',xdim,1,zdim,12);
RG.pressure_3D = repmat(permute(RG.pressure,[3 2 1]),xdim,ydim,1,12);

%% fit GMM from climatological mean temperature and salinity
% transform to normalized arrays
idx = ~isnan(RG.temp_clim) & ~isnan(RG.sal_clim);
[X_norm,C,S] = normalize([RG.temp_clim(idx) RG.sal_clim(idx) RG.lon_cos_3D(idx)...
    RG.latitude_3D(idx) RG.pressure_3D(idx)]);
RG = rmfield(RG,{'temp_clim' 'sal_clim' 'lon_cos_3D' 'latitude_3D' 'pressure_3D'});
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

%% assign grid cells and probabilities to clusters
% load GMM model
load(['Data/GMM_' num2str(num_clusters) '/model'],'gmm','C','S');
% for each month
cos_lon = repmat(cosd(RG.longitude-20),1,length(RG.latitude),length(RG.pressure));
lat = repmat(RG.latitude',length(RG.longitude),1,length(RG.pressure));
pres = repmat(permute(RG.pressure,[3 2 1]),length(RG.longitude),length(RG.latitude),1);
for m = 1:length(RG.time)
    % transform to normalized arrays
    temp = RG.temp(:,:,:,m);
    sal = RG.sal(:,:,:,m);
    idx = ~isnan(temp) & ~isnan(sal);
    X_norm = normalize([temp(idx) sal(idx) cos_lon(idx) lat(idx) pres(idx)],...
        'Center',C,'Scale',S);
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters
    GMM_monthly = nan(size(RG.temp,1),size(RG.temp,2),size(RG.temp,3));
    GMM_monthly(idx) = clusters;
    % save clusters
    if ~isfolder(['Data/GMM_' num2str(num_clusters)])
        mkdir(['Data/GMM_' num2str(num_clusters)]);
    end
    save(['Data/GMM_' num2str(num_clusters) '/m' num2str(m)],'GMM_monthly','-v7.3');
    clear GMM_monthly
    % fill 3D probabilities
    for c = 1:num_clusters
        GMM_monthly_probs = nan(size(RG.temp,1),size(RG.temp,2),size(RG.temp,3));
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
GMM.longitude = RG.longitude;
GMM.latitude = RG.latitude;
GMM.pressure = RG.pressure;
GMM.time = RG.time;
for m = 1:length(RG.time)
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
    GMM_probs.longitude = RG.longitude;
    GMM_probs.latitude = RG.latitude;
    GMM_probs.pressure = RG.pressure;
    GMM_probs.time = RG.time;
    for m = 1:length(RG.time)
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
clear RG

toc
% 1.8 hours for ten clusters and ten replicates on chinook? (9/8/23)
% 2.2 hours for twenty clusters and twenty replicates on Hercules (9/8/23)
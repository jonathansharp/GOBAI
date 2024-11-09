% gmm_clustering
%
% DESCRIPTION:
% This function uses gridded temperature and salinity with Gaussian Mixture
% Modelling to formulate global clusters of similar environmental
% variability within which to train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 6/18/2024

function gmm_clustering(fpath,base_grid,start_year,snap_date,clust_vars,num_clusters,numWorkers_predict)

%% process date
date_str = num2str(snap_date);

% check for existence of model
if exist(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model_' date_str '.mat'],'file') ~= 2

%% load climatological temperature and salinity data
[TS,timesteps] = load_climatological_TS_data(fpath,base_grid,start_year,date_str);

%% load chlorophyll climatological data
% if strcmp(base_grid,'RG')
%     load_chl();
% 
% elseif strcmp(base_grid,'RFROM')
%     fpath = 'Data/RFROM/RFROM_CHL_CMEMS_annual/';
%     chl = ncread(fpath,'ocean_temperature');
% 
% end

%% load basin mask file


%% fit GMM from climatological mean temperature and salinity
% Some replicates don't converge, investigate further...
tic
% transform to normalized arrays
idx = ~isnan(TS.temperature_cns) & ~isnan(TS.salinity_abs);
predictor_matrix = [];
for v = 1:length(clust_vars)
    predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
end
[X_norm,C,S] = normalize(predictor_matrix);
clear TS
% reduce inputs for model training to 100,000 random data points
idx_rand = randperm(length(X_norm),10000)';
X_norm = X_norm(idx_rand,:);
% fit GMM
gmm = fitgmdist(X_norm,num_clusters,...
    'CovarianceType','full',...
    'SharedCovariance',true,'Replicates',20);
% save GMM model
if ~isfolder(['Data/GMM_' base_grid '_' num2str(num_clusters)])
    mkdir(['Data/GMM_' base_grid '_' num2str(num_clusters)]);
end
save(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model_' date_str],...
    'gmm','num_clusters','C','S','-v7.3');
clear gmm C S
toc

%% assign grid cells and probabilities to clusters

% start timing cluster assignment
tic

% load GMM model
load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model_' ...
    date_str],'gmm','C','S');

% set up parallel pool
% tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

% for each month
for m = 1:timesteps

    % parse number of weeks in each month for RFROM
    if strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim([fpath '/Data/RFROM/']);
        % determine number of weeks in file
        nc_atts = ncinfo([fpath '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        max_w = nc_atts.Dimensions(3).Length;
    else
        max_w = 1;
    end

    for w = 1:max_w

        % load T/S grid
        TS = load_monthly_TS_data(fpath,base_grid,m,w,start_year,date_str);

        % transform to normalized arrays
        idx = ~isnan(TS.temperature_cns) & ~isnan(TS.salinity_abs);
        predictor_matrix = [];
        for v = 1:length(clust_vars)
            predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
        end
        X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
        % assign data points to clusters
        assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,TS.xdim,...
            TS.ydim,TS.zdim,m,1,base_grid);

    end

end

% end parallel session
delete(gcp('nocreate'));

toc

% display information
disp([num2str(num_clusters) ' clusters formed using ' base_grid ' grid']);

else

% display information
disp(['already used ' base_grid ' grid to form ' num2str(num_clusters) ...
    ' clusters for ' date_str(5:6) '/' date_str(1:4)]);

end

end

%% embedded function to assign points to GMM clusters
function assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,xdim,ydim,zdim,m,w,base_grid)
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters (highest probability cluster)
    GMM_clusters = nan(xdim,ydim,zdim);
    GMM_clusters(idx) = clusters;
    % save clusters
    folder_name = [pwd '/Data/GMM_' base_grid '_' num2str(num_clusters)];
    if ~isfolder(folder_name); mkdir(folder_name); end
    save([folder_name '/m' num2str(m) '_w' num2str(w)],'GMM_clusters','-v7.3');
    clear GMM_clusters
    % fill 3D probabilities (for each cluster)
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

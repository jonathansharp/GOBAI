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

function gmm_clustering(base_grid,clust_vars,num_clusters,numWorkers_predict)

% check for existence of model
if exist(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model.mat'],'file') ~= 2

%% load climatological data
if strcmp(base_grid,'RG')
    % load practical salinity and in situ temperature
    [TS,timesteps] = load_RG_dim([pwd '/Data/RG_CLIM/']);
    TS = load_RG_clim(TS,[pwd '/Data/RG_CLIM/']);
    % calculate absolute salinity and conservative temperature
    TS = replicate_RG_dim(TS,12);
    idx = ~isnan(TS.salinity) & ~isnan(TS.temperature);
    TS.salinity_abs = nan(size(idx));
    TS.salinity_abs(idx) = gsw_SA_from_SP(TS.salinity(idx),TS.pressure(idx),...
        convert_lon(TS.longitude(idx)),TS.latitude(idx));
    TS.temperature_cns = nan(size(idx));
    TS.temperature_cns(idx) = ...
        gsw_CT_from_t(TS.salinity_abs(idx),TS.temperature(idx),TS.pressure(idx));
    TS = rmfield(TS,{'temperature' 'salinity'});
    % load chlorophyll

elseif strcmp(base_grid,'RFROM')
    % load absolute salinity and conservative temperature
    [TS,timesteps] = load_RFROM_dim([pwd '/Data/RFROM/']);
    TS = load_RFROM_clim(TS,[pwd '/Data/RFROM/']);
    TS = replicate_RFROM_dim(TS,1);
    % load chlorophyll 
end

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
save(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model'],...
    'gmm','num_clusters','C','S','-v7.3');
clear gmm C S
toc

%% assign grid cells and probabilities to clusters

% start timing cluster assignment
tic

% load GMM model
load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model'],'gmm','C','S');

% set up parallel pool
%p = setup_pool(numWorkers_predict);

% for each month
for m = 1%1:timesteps
    % load T/S grid
    if strcmp(base_grid,'RG')
        % load dimensions
        TS = load_RG_dim([pwd '/Data/RG_CLIM/']);
        TS = replicate_RG_dim(TS,1);
        % get RG practical salinity and in situ temperature
        TS.temperature = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],...
            'Temperature',[1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Sal.nc'],...
            'Salinity',[1 1 1 m],[Inf Inf Inf 1]);
        % convert to absolute salinity and conservative temperature
        idx = ~isnan(TS.salinity) & ~isnan(TS.temperature);
        TS.salinity_abs = nan(size(idx));
        TS.salinity_abs(idx) = gsw_SA_from_SP(TS.salinity(idx),TS.pressure(idx),...
            convert_lon(TS.longitude(idx)),TS.latitude(idx));
        TS.temperature_cns = nan(size(idx));
        TS.temperature_cns(idx) = ...
            gsw_CT_from_t(TS.salinity_abs(idx),TS.temperature(idx),TS.pressure(idx));
        % transform to normalized arrays
        idx = ~isnan(TS.temperature) & ~isnan(TS.salinity);
        predictor_matrix = [];
        for v = 1:length(clust_vars)
            predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
        end
        X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
        % assign data points to clusters
        assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,TS.xdim,...
            TS.ydim,TS.zdim,m,1,base_grid);
    else
        % load dimensions
        TS = load_RFROM_dim([pwd '/Data/RFROM/']);
        TS = replicate_RFROM_dim(TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % get RFROM absolute salinity and conservative temperature
            TS.temperature_cns = ncread([pwd '/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity_abs = ncread([pwd '/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % transform to normalized arrays
            idx = ~isnan(TS.temperature_cns) & ~isnan(TS.salinity_abs);
            predictor_matrix = [];
            for v = 1:length(clust_vars)
                predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
            end
            X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
            % assign data points to clusters
            assign_to_gmm_clusters(gmm,num_clusters,idx,X_norm,...
                TS.xdim,TS.ydim,TS.zdim,m,w,base_grid);
        end
    end
end

% end parallel session
delete(gcp('nocreate'));

toc

% display information
disp([num2str(num_clusters) ' clusters formed using ' base_grid ' grid']);

else

% display information
disp(['already used ' base_grid ' grid to form ' num2str(num_clusters) ' clusters']);

end

end

%% embedded function
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

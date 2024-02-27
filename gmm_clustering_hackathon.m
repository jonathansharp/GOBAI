% gmm_clustering_hackathon
%
% DESCRIPTION:
% This function uses gridded temperature and salinity with Gaussian Mixture
% Modelling to formulate global clusters of similar environmental
% variability within which to train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 12/1/2023

%% initiate profile
profile on

%% load configuration parameters
gobai_o2_initiate;
load_standard_config_files;
load('Config/base_config_RFROM.mat'); % base grid
load('Config/predict_years_config_04.mat'); % only for 2004
dir_base = create_dir_base('FFNN',{base_grid;num_clusters;file_date;...
    float_file_ext;train_ratio;val_ratio;test_ratio}); % directory name base
fpath = '/raid'; % for RFROM
% fpath = pwd; % for RG

%% load temperature and salinity data
if strcmp(base_grid,'RG')
    [TS,timesteps] = load_RG_dim([pwd '/Data/RG_CLIM/']);
    TS = load_RG_clim(TS,[pwd '/Data/RG_CLIM/']);
    TS = replicate_RG_dim(TS,12);
elseif strcmp(base_grid,'RFROM')
    [TS,timesteps] = load_RFROM_dim('/raid/Data/RFROM/');
    TS = load_RFROM_clim(TS,'/raid/Data/RFROM/');
    TS = replicate_RFROM_dim(TS,1);
end

%% fit GMM from climatological mean temperature and salinity
tic
% transform to normalized arrays
idx = ~isnan(TS.temperature) & ~isnan(TS.salinity);
predictor_matrix = [];
for v = 1:length(clust_vars)
    predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
end
[X_norm,C,S] = normalize(predictor_matrix);
clear TS
% reduce inputs from 30,000,000 to 100,000 random data points %%%%%%%%%%%%%
idx_rand = randperm(length(X_norm),100000)';
% fit GMM
gmm = fitgmdist(X_norm(idx_rand,:),num_clusters,...
    'CovarianceType','full',...
    'SharedCovariance',true,'Replicates',10);
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
%tic; parpool(12); fprintf('Pool initiation:'); toc;

% for each month (1)
for m = 1
    % load T/S grid
    if strcmp(base_grid,'RG')
        % load dimensions
        TS = load_RG_dim([pwd '/Data/RG_CLIM/']);
        TS = replicate_RG_dim(TS,1);
        % get RG T/S
        TS.temperature = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Temp.nc'],...
            'Temperature',[1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([pwd '/Data/RG_CLIM/RG_Climatology_Sal.nc'],...
            'Salinity',[1 1 1 m],[Inf Inf Inf 1]);
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
        TS = load_RFROM_dim('/raid/Data/RFROM/');
        TS = replicate_RFROM_dim(TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo(['/raid/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        % just week 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for w = 1%:nc_atts.Dimensions(3).Length
            % get RFROM T and S
            TS.temperature = ncread(['/raid/Data/RFROM/RFROM_TEMP_v0.1/RFROM_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity = ncread(['/raid/Data/RFROM/RFROM_SAL_v0.1/RFROM_SAL_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % transform to normalized arrays
            idx = ~isnan(TS.temperature) & ~isnan(TS.salinity);
            % just take top layer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idx(:,:,2:end) = [];
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

%% end and save profile
p=profile('info');
profsave(p,'profiles/gmm_clustering_hackathon')
profile off

%% embedded function
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

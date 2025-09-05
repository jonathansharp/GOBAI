% gmm_clustering
%
% DESCRIPTION:
% This function uses gridded temperature and salinity with Gaussian Mixture
% Modelling to formulate global clusters of similar environmental
% variability within which to train machine learning models.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/31/2025

function gmm_clustering_new(param_props,fpaths,base_grid,start_year,end_year,snap_date,...
    float_file_ext,clust_vars,num_clusters,numWorkers_cluster,flt,gld,ctd,varargin)

%% define dataset extensions
if flt == 1; float_ext = 'f'; else float_ext = ''; end
if gld == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

%% process optional input arguments
% pre-allocate
rlz = NaN;
% process inputs
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'rlz')
        rlz = varargin{i+1};
    end
end

%% process date
date_str = num2str(snap_date);
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');

%% define GMM model name
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
        num2str(num_clusters)];
    if ~isfolder(gmm_folder_name); mkdir(gmm_folder_name); end
    gmm_model_name = [gmm_folder_name '/model_' float_ext ...
        glodap_ext ctd_ext '_' date_str];
else
    gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
        base_grid '_' num2str(num_clusters)];
    if ~isfolder(gmm_folder_name); mkdir(gmm_folder_name); end
    gmm_model_name = [gmm_folder_name '/model_' float_ext ...
        glodap_ext ctd_ext '_' date_str];
end

% %% check for existence of model
% if ~isfile([gmm_model_name '.mat'])

    %% load basin mask file
    
    
    %% load data
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        load([param_props.dir_name '/Data/processed_all_' param_props.file_name ...
            '_data_' float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],...
             'all_data');
    else
        load([param_props.dir_name '/Data/' base_grid '_' param_props.file_name '_data_' ...
            float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],'all_data');
    end
    
    %% fit GMM from data points themselves
    tic
    % transform to normalized arrays
    idx = false(length(all_data.temperature),1);
    for v = 1:length(clust_vars)
        idx = idx | ~isnan(all_data.(clust_vars{v}));
    end
    predictor_matrix = [];
    for v = 1:length(clust_vars)
        predictor_matrix = [predictor_matrix all_data.(clust_vars{v})];
    end
    basins = find_basin_paige();
    for b = 1
    [X_norm,C,S] = normalize(predictor_matrix);
    % reduce inputs for model training to 100,000 random data points
    idx_rand = randperm(length(X_norm),100000)';
    X_norm = X_norm(idx_rand,:);
    % fit GMM
    options = statset('MaxIter',1000); % increase max iterations to ensure convergence
    gmm = fitgmdist(X_norm,num_clusters,...
        'Options',options,...
        'CovarianceType','diagonal',...
        'SharedCovariance',true,'Replicates',100);
    % save GMM model
    save(gmm_model_name,'gmm','num_clusters','C','S','-v7.3');

    % plot clusters by temperature and salinity
    % load([param_props.dir_name '/Data/GMM_' num2str(num_clusters) '/model_' date_str],...
    %      'gmm','num_clusters','C','S');
    % clusters = cluster(gmm,X_norm);
    % figure; hold on;
    % predictor_matrix_reduced = predictor_matrix(idx_rand,:);
    % h = gscatter(predictor_matrix_reduced(:,2),...
    %     predictor_matrix_reduced(:,1),clusters,flipud(jet(num_clusters)));
    % colormap();
    % xlabel('Salinity');
    % ylabel('Temperature');
    % title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
    % hold off;
    % export('');
    % close


    % clean up
    clear gmm C S
    toc
    
    %% fit GMM from climatological mean temperature and salinity
    
    % load climatological temperature and salinity data
    % TS = load_climatological_TS_data(fpaths.temp_path,base_grid,start_year,date_str);
    % 
    % % Some replicates don't converge, investigate further...
    % tic
    % % transform to normalized arrays
    % idx = ~isnan(TS.temperature_cns) & ~isnan(TS.salinity_abs);
    % predictor_matrix = [];
    % for v = 1:length(clust_vars)
    %     predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
    % end
    % [X_norm,C,S] = normalize(predictor_matrix);
    % clear TS
    % % reduce inputs for model training to 100,000 random data points
    % idx_rand = randperm(length(X_norm),100000)';
    % X_norm = X_norm(idx_rand,:);
    % % fit GMM
    % options = statset('MaxIter',1000); % increase max iterations to ensure convergence
    % gmm = fitgmdist(X_norm,num_clusters,...
    %     'Options',options,...
    %     'CovarianceType','diagonal',...
    %     'SharedCovariance',true,'Replicates',10);
    % % save GMM model
    % if ~isfolder(['Data/GMM_' base_grid '_' num2str(num_clusters)])
    %     mkdir(['Data/GMM_' base_grid '_' num2str(num_clusters)]);
    % end
    % save(['Data/GMM_' base_grid '_' num2str(num_clusters) '/model_' date_str],...
    %     'gmm','num_clusters','C','S','-v7.3');
    % clear gmm C S
    % toc
% 
% else
% 
%     % display information
%     disp(['already trained GMM for ' num2str(num_clusters) ...
%         ' clusters for ' date_str(5:6) '/' date_str(1:4)]);
% 
% end

%% assign grid cells and probabilities to clusters
folder_name = [fpaths.param_path 'GMM_' base_grid '_' num2str(num_clusters)];
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    folder_name_temp = [fpaths.param_path_temp 'GMM_' base_grid '_' num2str(num_clusters)];
else
    folder_name_temp = [fpaths.model_path base_grid '/GMM_' base_grid '_' num2str(num_clusters)];
end
% determine length of cluster file if it exists
if exist([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'file') == 2
    inf = ncinfo([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc']);
    for n = 1:length(inf.Dimensions)
        if strcmp(inf.Dimensions(n).Name,'time')
            time_idx = n;
        end
    end
end
% determine expected length
% year = str2double(date_str(1:4));
% month = str2double(date_str(5:6));
% length_expt = (year-start_year)*12 + month;
length_expt = (end_year-start_year)*12 + 12;

%% check for existence of cluster file and length of cluster grids
% if ~isfile([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc']) || ...
%         inf.Dimensions(time_idx).Length ~= length_expt

    % delete file if it exists
    if exist([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'],'file') == 2
        delete([folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc']);
    end

    % start timing cluster assignment
    tic
    
    %% create netCDF file that will be end result
    if strcmp(base_grid,'RG')
       [TS,months,weeks,timesteps] = load_RG_dim(fpaths.temp_path);
       % create file
       create_nc_files(TS,num_clusters,base_grid,TS.xdim,TS.ydim,TS.zdim,folder_name);
    elseif strcmp(base_grid,'RFROM')
        [TS,months,weeks,timesteps] = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
        % create file
        create_nc_files(TS,0,base_grid,TS.xdim,TS.ydim,TS.zdim,folder_name); % don't write cluster probs for RFROM
    else
        % define paths
        path2 = ['_Omon_' base_grid '_'];
        path3 = ['_' rlz '_gr'];
        % define filepaths
        nc_filepath_abs_sal = [fpaths.model_path base_grid ...
            '/combined/regridded/abs_sal' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        % load dimensions
        [TS,months,weeks,timesteps] = load_model_dim(nc_filepath_abs_sal);
        % create file
        create_nc_files(TS,num_clusters,base_grid,TS.xdim,TS.ydim,TS.zdim,folder_name);
    end
    
    % load GMM model
    load(gmm_model_name,'gmm','C','S');
    
    % set up parallel pool
    tic; parpool(numWorkers_cluster); fprintf('Pool initiation: '); toc;

    %% assign clusters for each timestep
    parfor t = 1:timesteps

        % load T/S grid
        TS = load_monthly_TS_data(fpaths,base_grid,...
            months(t),weeks(t),start_year,end_year,date_str);

        % transform to normalized arrays
        idx = ~isnan(TS.temperature_cns) & ~isnan(TS.salinity_abs);
        predictor_matrix = [];
        for v = 1:length(clust_vars)
            predictor_matrix = [predictor_matrix TS.(clust_vars{v})(idx)];
        end
        X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);

        % assign data points to clusters
        assign_to_gmm_clusters(TS,base_grid,gmm,num_clusters,idx,X_norm,t,folder_name_temp);

    end

    % end parallel session
    delete(gcp('nocreate'));

    %% concatenate cluster information in main file
    for t = 1:timesteps
        % define file names
        filename = [folder_name '/clusters_' float_ext glodap_ext ctd_ext '.nc'];
        filename_temp = [folder_name_temp '/clust_' num2str(t) '.nc'];
        % read information from temporary file and write it to main file
        time = ncread(filename_temp,'time'); % read
        ncwrite(filename,'time',time,t); % write
        GMM_clusters = ncread(filename_temp,'GMM_clusters'); % read
        ncwrite(filename,'clusters',GMM_clusters,[1 1 1 t]); % write
        if strcmp(base_grid,'RFROM')
            % do nothing
        else
            for c = 1:num_clusters
                GMM_cluster_probs = ncread(filename_temp,...
                    ['GMM_cluster_probs_' num2str(c)]); % read
                ncwrite(filename,['cluster_probs_c' num2str(c)],...
                    GMM_cluster_probs,[1 1 1 t]); % write
            end
        end
        % delete temporary file
        delete(filename_temp);
    end
    
    % display information
    toc
    disp([num2str(num_clusters) ' clusters formed using ' base_grid ' grid']);

% else
% 
%     % display information
%     disp(['already used ' base_grid ' grid to form ' num2str(num_clusters) ...
%         ' clusters for ' date_str(5:6) '/' date_str(1:4)]);
% 
% end

end

%% embedded function to assign points to GMM clusters
function assign_to_gmm_clusters(TS,base_grid,gmm,num_clusters,idx,X_norm,t,folder_name_temp)
    % assign to clusters and obtain probabilities
    [clusters,~,p] = cluster(gmm,X_norm);
    % fill 3D clusters (highest probability cluster)
    GMM_clusters = nan(TS.xdim,TS.ydim,TS.zdim);
    GMM_clusters(idx) = uint8(clusters);
    % save cluster properties in temporary files
    filename = [folder_name_temp '/clust_' num2str(t) '.nc'];
    if isfile(filename); delete(filename); end
    nccreate(filename,'time','Dimensions',{'time' 1});
    ncwrite(filename,'time',TS.Time(t));
    nccreate(filename,'GMM_clusters','Dimensions',{'lon' size(GMM_clusters,1) ...
        'lat' size(GMM_clusters,2) 'pres' size(GMM_clusters,3)});
    ncwrite(filename,'GMM_clusters',GMM_clusters);
    if strcmp(base_grid,'RFROM')
        % do nothing
    else
        % fill 3D probabilities (for each cluster)
        for c = 1:num_clusters
            GMM_cluster_probs = nan(TS.xdim,TS.ydim,TS.zdim);
            GMM_cluster_probs(idx) = int16(p(:,c)*10000);
            nccreate(filename,['GMM_cluster_probs_' num2str(c)],...
                'Dimensions',{'lon' size(GMM_cluster_probs,1) 'lat' ...
                size(GMM_cluster_probs,2) 'pres' size(GMM_cluster_probs,3)});
            ncwrite(filename,['GMM_cluster_probs_' num2str(c)],...
                GMM_cluster_probs);
    end
    end
end

% for creating netCDF file
function create_nc_files(TS,num_clusters,base_grid,xdim,ydim,zdim,folder_name)

    % define file name and create file
    if ~isfolder(folder_name); mkdir(folder_name); end
    filename = [folder_name '/clusters.nc'];

    % longitude
    nccreate(filename,'lon','Dimensions',{'lon',xdim},...
        'DataType','single','FillValue',NaN,'Format','netcdf4');
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
    
    % pressure (or depth)
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        nccreate(filename,'pres','Dimensions',{'pres',zdim},...
            'DataType','single','FillValue',NaN);
        ncwrite(filename,'pres',TS.Pressure);
        ncwriteatt(filename,'pres','units','decibars');
        ncwriteatt(filename,'pres','axis','Z');
        ncwriteatt(filename,'pres','long_name','pressure');
        ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');
    else
        nccreate(filename,'depth','Dimensions',{'depth',zdim},...
        'DataType','single','FillValue',NaN);
        ncwrite(filename,'depth',TS.Depth);
        ncwriteatt(filename,'depth','units','meters');
        ncwriteatt(filename,'depth','axis','Z');
        ncwriteatt(filename,'depth','long_name','depth');
        ncwriteatt(filename,'depth','_CoordinateAxisType','Depth');
    end
    
    % time
    nccreate(filename,'time','Dimensions',{'time',Inf},...
        'DataType','single','FillValue',NaN);
    ncwriteatt(filename,'time','units','days since 0000-01-01');
    ncwriteatt(filename,'time','axis','T');
    ncwriteatt(filename,'time','long_name','time');
    ncwriteatt(filename,'time','_CoordinateAxisType','Time');

    % clusters
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        nccreate(filename,'clusters','Dimensions',{'lon',xdim,'lat',ydim,...
            'pres',zdim,'time',Inf},'DataType','uint8','FillValue',NaN);
        for c = 1:num_clusters
            nccreate(filename,['cluster_probs_c' num2str(c)],'Dimensions',{'lon',xdim,'lat',ydim,...
                'pres',zdim,'time',Inf},'DataType','int16','FillValue',NaN);
        end
    else
        nccreate(filename,'clusters','Dimensions',{'lon',xdim,'lat',ydim,...
            'depth',zdim,'time',Inf},'DataType','uint8','FillValue',NaN);
        for c = 1:num_clusters
            nccreate(filename,['cluster_probs_c' num2str(c)],'Dimensions',{'lon',xdim,'lat',ydim,...
                'depth',zdim,'time',Inf},'DataType','int16','FillValue',NaN);
        end
    end

end

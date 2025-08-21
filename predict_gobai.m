% predict_gobai
%
% DESCRIPTION:
% This function uses trained models to predict
% on a grid.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 7/31/2025

function predict_gobai(alg_type,param_props,fpaths,base_grid,file_date,...
    float_file_ext,num_clusters,variables,thresh,numWorkers_predict,...
    clust_vars,start_year,end_year,snap_date,flt,gld,ctd,varargin)

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

%% process necessary input arguments for model parameters
% pre-allocate
train_ratio = NaN;
val_ratio = NaN;
test_ratio = NaN;
numtrees = NaN;
minLeafSize = NaN;
numstumps = NaN;
numbins = NaN;
% process based on algorithm type
if strcmp(alg_type,'FFNN')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'train_ratio')
            train_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'val_ratio')
            val_ratio = varargin{i+1};
        elseif strcmpi(varargin{i}, 'test_ratio')
            test_ratio = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'RFR')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numtrees')
            numtrees = varargin{i+1};
        elseif strcmpi(varargin{i}, 'minLeafSize')
            minLeafSize = varargin{i+1};
        end
    end
elseif strcmp(alg_type,'GBM')
    for i = 1:2:length(varargin)-1
        if strcmpi(varargin{i}, 'numstumps')
            numstumps = varargin{i+1};
        elseif strcmpi(varargin{i}, 'numbins')
            numbins = varargin{i+1};
        end
    end
else
    disp('"alg_type" must be "FFNN", "RFR", or "GBM"')
end

%% process date
date_str = num2str(snap_date);

%% directory base
if strcmp(alg_type,'FFNN')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
elseif strcmp(alg_type,'RFR')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
elseif strcmp(alg_type,'GBM')
    dir_base = create_dir_base(alg_type,{num_clusters;file_date;...
        float_file_ext;numstumps;numbins});
end

%% create directory and file names
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    alg_dir = [param_props.dir_name '/Models/' dir_base];
else
    alg_dir = [param_props.dir_name '/Models/' base_grid '/' dir_base];
end
alg_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    alg_fnames(c) = ...
        {[alg_type '_' param_props.file_name '_C' num2str(c) '_' float_ext glodap_ext ctd_ext]};
end
if strcmp(alg_type,'FFNN')
    path_ext = ['GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/train' num2str(100*train_ratio) ...
        '_val' num2str(100*val_ratio) '_test' num2str(100*test_ratio) ...
        '/' float_ext glodap_ext ctd_ext '/'];
elseif strcmp(alg_type,'RFR')
    path_ext = ['GOBAI/' base_grid '/RFR/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numtrees) '_lf' ...
        num2str(minLeafSize) '/' float_ext glodap_ext ctd_ext '/'];
elseif strcmp(alg_type,'GBM')
    path_ext = ['GOBAI/' base_grid '/GBM/c' num2str(num_clusters) ...
        '_' file_date float_file_ext '/tr' num2str(numstumps) ...
        '_bin' num2str(numbins) '/' float_ext glodap_ext ctd_ext '/'];
end
gobai_alg_dir_temp = [fpaths.param_path_temp path_ext];
if ~isfolder(gobai_alg_dir_temp); mkdir(gobai_alg_dir_temp); end
gobai_alg_dir = [fpaths.param_path path_ext];
if ~isfolder(gobai_alg_dir); mkdir(gobai_alg_dir); end

%% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% load dimensions
if strcmp(base_grid,'RG')
    TS = load_RG_dim(fpaths.temp_path);
elseif strcmp(base_grid,'RFROM')
    TS = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = ['_' rlz '_gr'];
    % define filepaths
    nc_filepath_abs_sal = [fpaths.model_path base_grid '/combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
end

%% set up parallel pool
tic; parpool(numWorkers_predict); fprintf('Pool initiation: '); toc;

%% start timing predictions
tStart = tic;

%% compute and save estimates for each month
parfor m = 1:length(TS.months)
    if strcmp(base_grid,'RG')
        % counter 
        cnt = m;
        % load dimensions
        TS = load_RG_dim(fpaths.temp_path);
        TS = replicate_dims(base_grid,TS,1);
        TS.longitude = convert_lon(TS.longitude);
        % get RG T and S
        TS.temperature = ncread([fpaths.temp_path 'RG_Climatology_Temp.nc'],...
            'Temperature',[1 1 1 m],[Inf Inf Inf 1]);
        TS.salinity = ncread([fpaths.sal_path 'RG_Climatology_Sal.nc'],...
            'Salinity',[1 1 1 m],[Inf Inf Inf 1]);
        % covert RG T and S to conservative temperature and absolute salinity
        pres_3d = repmat(permute(TS.Pressure,[3 2 1]),length(TS.Longitude),length(TS.Latitude),1);
        lon_3d = repmat(TS.Longitude,1,length(TS.Latitude),length(TS.Pressure));
        lat_3d = repmat(TS.Latitude',length(TS.Longitude),1,length(TS.Pressure));
        TS.salinity_abs = gsw_SA_from_SP(TS.salinity,pres_3d,convert_lon(lon_3d),lat_3d);
        TS.temperature_cns = gsw_CT_from_t(TS.salinity_abs,TS.temperature,pres_3d);
        % get time variables for just this timestep
        TS.Time = ncread([fpaths.temp_path 'RG_Climatology_Temp.nc'],'Time',m,1);
        date_temp = datevec(datenum(2004,1,1+double(TS.Time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply model
        apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
            base_grid,m,1,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_alg_dir_temp,param_props,fpaths.param_path,...
            date_str,clust_vars,float_ext,glodap_ext,ctd_ext);
    elseif strcmp(base_grid,'RFROM')
        % load dimensions
        TS = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
        TS = replicate_dims(base_grid,TS,1);
        % determine number of weeks in file
        nc_atts = ncinfo([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
            num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc']);
        for w = 1:nc_atts.Dimensions(3).Length
            % counter
            cnt = TS.cnt{m}(w);
            % get RFROM T and S
            TS.temperature_cns = ncread([fpaths.temp_path 'RFROM_TEMP_v2.2/RFROMV22_TEMP_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                    'ocean_temperature',[1 1 1 w],[Inf Inf Inf 1]);
            TS.salinity_abs = ncread([fpaths.sal_path 'RFROM_SAL_v2.2/RFROMV22_SAL_STABLE_' ...
                num2str(TS.years(m)) '_' sprintf('%02d',TS.months(m)) '.nc'],...
                    'ocean_salinity',[1 1 1 w],[Inf Inf Inf 1]);
            % load bgc variables if applicable
            if any(strcmp(variables,'o2'))
                TS.o2 = ncread(['/fast4/o2/GOBAI/' base_grid '/' dir_base '/gobai-o2.nc'],...
                        'o2',[1 1 1 cnt],[Inf Inf Inf 1]);
            end
            if any(strcmp(variables,'no3'))
                TS.no3 = ncread(['/fast5/no3/GOBAI/' base_grid '/' dir_base '/gobai-no3.nc'],...
                        'no3',[1 1 1 cnt],[Inf Inf Inf 1]);
            end
            % get time variables for just this timestep
            date_temp = datevec(datenum(1950,0,0)+TS.Time(cnt));
            date_temp0 = date_temp;
            date_temp0(:,2:3) = 1; % Jan. 1 of the year
            TS.year = date_temp(:,1);
            TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
            % transform day
            TS.day_sin = sin((2.*pi.*TS.day)/365.25);
            TS.day_cos = cos((2.*pi.*TS.day)/365.25);
            % apply model
            apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
                base_grid,m,w,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
                thresh,gobai_alg_dir_temp,param_props,fpaths.param_path,...
                date_str,clust_vars,float_ext,glodap_ext,ctd_ext);
        end
    else
        % counter 
        cnt = m;
        % define paths
        path2 = ['_Omon_' base_grid '_'];
        path3 = ['_'  rlz '_gr'];
        % define filepaths
        nc_filepath_abs_sal = [fpaths.model_path base_grid '/combined/regridded/abs_sal' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        nc_filepath_cns_tmp = [fpaths.model_path base_grid '/combined/regridded/cns_tmp' path2 ...
            'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
        % load dimensions
        TS = load_model_dim(nc_filepath_abs_sal);
        TS = replicate_dims(base_grid,TS,1);
        % get practical salinity and in situ temperature from cmip model
        TS.salinity_abs = ncread(nc_filepath_abs_sal,'abs_sal',[1 1 1 m],[Inf Inf Inf 1]);
        TS.temperature_cns = ncread(nc_filepath_cns_tmp,'cns_tmp',[1 1 1 m],[Inf Inf Inf 1]);
        % get time variables for just this timestep
        TS.Time = ncread(nc_filepath_abs_sal,'time',m,1);
        date_temp = datevec(datenum(0,0,double(TS.Time)));
        date_temp0 = date_temp;
        date_temp0(:,2:3) = 1; % Jan. 1 of each year
        TS.year = date_temp(:,1);
        TS.day = datenum(date_temp) - datenum(date_temp0) + 1;
        % transform day
        TS.day_sin = sin((2.*pi.*TS.day)/365.25);
        TS.day_cos = cos((2.*pi.*TS.day)/365.25);
        % apply model
        apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
            base_grid,m,1,cnt,TS.xdim,TS.ydim,TS.zdim,variables_TS,...
            thresh,gobai_alg_dir_temp,param_props,fpaths.param_path,...
            date_str,clust_vars,float_ext,glodap_ext,ctd_ext);

    end
end

% end parallel session
delete(gcp('nocreate'));

%% create netCDF file that will be end result
if strcmp(base_grid,'RG')
    TS = load_RG_dim(fpaths.temp_path);
   % create file
   create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param_props);
elseif strcmp(base_grid,'RFROM')
    TS = load_RFROM_dim(fpaths.temp_path,start_year,end_year);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param_props);
else
    % define paths
    path2 = ['_Omon_' base_grid '_'];
    path3 = ['_' rlz '_gr'];
    % define filepaths
    nc_filepath_abs_sal = [fpaths.model_path base_grid '/combined/regridded/abs_sal' path2 ...
        'combined' path3 '_' num2str(start_year) '01-' date_str '.nc'];
    % load dimensions
    TS = load_model_dim(nc_filepath_abs_sal);
    % create file
    create_nc_file(TS,base_grid,TS.xdim,TS.ydim,TS.zdim,gobai_alg_dir,param_props);
end

%% concatenate gobai in main file
files = dir([gobai_alg_dir_temp '/gobai-' param_props.file_name '-*.nc']); % count files in folder
idx_rem = [];
for fln = 1:length(files)
    % index to remove filenames with uncertinaty
    if contains(files(fln).name,'uncer')
        idx_rem = [idx_rem;fln];
    end
end
files(idx_rem) = [];
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];
for cnt = 1:length(files)
    % define file name
    filename_temp = [gobai_alg_dir_temp 'gobai-' ...
        param_props.file_name '-' num2str(cnt) '.nc'];
    % read information from temporary file and write it to main file
    time = ncread(filename_temp,'time'); % read
    ncwrite(filename,'time',time,cnt); % write
    gobai_3d = ncread(filename_temp,param_props.file_name); % read
    ncwrite(filename,param_props.file_name,gobai_3d,[1 1 1 cnt]); % write
    % delete temporary file
    delete(filename_temp);
end

%% create monthly,one degree file
% for x = 1:length(TS.Longitude)
%     for y = 1:length(TS.Latitude)
%         for z = 1:length(TS.Pressure)
%             var = squeeze(ncread(filename,param_props.file_name,[x y z 1],[1 1 1 Inf]));
%             if any(~isnan(var))
% 
%             end
%         end
%     end
% end

%% stop timing predictions
fprintf([alg_type ' Prediction (' num2str(start_year) ' to ' date_str(1:4) '): ']);

tElapsed = toc(tStart);
disp(['Elapsed time is ' num2str(tElapsed/60) ' minutes.'])

end

%% for processing 3D grids and applying trained models to them
function apply_model(alg_type,TS,num_clusters,alg_dir,alg_fnames,...
    base_grid,m,w,cnt,xdim,ydim,zdim,variables_TS,thresh,gobai_alg_dir_temp,...
    param_props,param_path,date_str,clust_vars,float_ext,glodap_ext,ctd_ext)

    % convert to arrays
    TS_index = ~isnan(TS.temperature_cns);
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
        if isscalar(TS.(vars{v}))
            TS.([vars{v} '_array']) = repmat(TS.(vars{v}),size(TS.temperature_cns_array));
            TS = rmfield(TS,vars{v});
        end
    end

    % calculate absolute salinity, conservative temperature, potential density
    TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

    % pre-allocate
    gobai_matrix = single(nan(length(TS.temperature_cns_array),num_clusters));
    probs_matrix = single(nan(length(TS.temperature_cns_array),num_clusters));

    % define GMM model name
    if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
        gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
            num2str(num_clusters)];
        gmm_model_name = [gmm_folder_name '/model_' float_ext ...
            glodap_ext ctd_ext '_' date_str];
    else
        gmm_folder_name = [param_props.dir_name '/Data/GMM_' ...
            base_grid '_' num2str(num_clusters)];
        gmm_model_name = [gmm_folder_name '/model_' float_ext ...
            glodap_ext ctd_ext '_' date_str];
    end
    % load GMM model
    load(gmm_model_name,'gmm','C','S');
    % transform to normalized arrays
    predictor_matrix = [];
    for v = 1:length(clust_vars)
        predictor_matrix = [predictor_matrix TS.([clust_vars{v} '_array'])];
    end
    X_norm = normalize(predictor_matrix,'Center',C,'Scale',S);
    % apply model to assign to clusters and obtain probabilities
    [~,~,p] = cluster(gmm,X_norm);
    clear predictor_matrix gmm X_norm

    % apply models for each cluster
    for c = 1:num_clusters

      % check for existence of model
      if isfile([alg_dir '/' alg_fnames{c} '.mat'])

        % load GMM cluster probabilities for this cluster and month, and convert to array
        %if strcmp(base_grid,'RFROM') % use 'p' for RFROM basegrid
            GMM_probs = nan(size(TS_index));
            GMM_probs(TS_index) = int16(p(:,c)*10000);
        %else % load probabilities for other basegrids %%%%%%% SOME ISSUE HERE
        %    GMM_probs = ncread([folder_name '/clusters.nc'],['cluster_probs_c' num2str(c)],...
        %        [1 1 1 cnt],[Inf Inf Inf 1]);
        %end
        probabilities_array = GMM_probs(TS_index)./10000; % convert to array
        probabilities_array(probabilities_array < thresh) = NaN; % remove probabilities below thresh
        probs_matrix(:,c) = probabilities_array; % add to probability matrix

        % load model for this cluster
        alg_struct = load([alg_dir '/' alg_fnames{c}],alg_type);
    
        % % remove year for cluster 13 for RG grid
        % variables_TS_tmp = variables_TS;
        % if c == 13
        %     variables_TS_tmp(strcmp(variables_TS_tmp,'year_array'))=[];
        % end

        % predict data for each cluster
        if strcmp(alg_type,'FFNN')
            gobai_matrix(:,c) = ...
                run_FFNN(alg_struct.FFNN,TS,probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS_tmp,thresh);
        elseif strcmp(alg_type,'RFR')
            gobai_matrix(:,c) = ...
                run_RFR(alg_struct.RFR,TS,probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS_tmp,thresh);
        elseif strcmp(alg_type,'GBM')
            gobai_matrix(:,c) = ...
                run_GBM(alg_struct.GBM,TS,probabilities_array,...
                true(size(TS.temperature_cns_array)),variables_TS_tmp,thresh);
        end

      end

    end

    % calculate weighted average over each cluster using probabilities
    gobai_array = ...
        double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
        sum(probs_matrix,2,'omitnan'));

    % zero trap
    gobai_array(gobai_array<0) = 0;
    
    % convert back to 3D grid
    gobai_3d = nan(xdim,ydim,zdim);
    gobai_3d(TS_index) = gobai_array;

    % Write output in temporary files
    filename = [gobai_alg_dir_temp 'gobai-' param_props.file_name '-' num2str(cnt) '.nc'];
    if exist(filename,'file')==2; delete(filename); end
    nccreate(filename,'time','Dimensions',{'time' 1});
    if strcmp(base_grid,'RFROM')
        ncwrite(filename,'time',TS.Time(cnt));
    else
        ncwrite(filename,'time',datenum(TS.years(cnt),TS.months(cnt),15)-datenum(1950,1,1));
    end
    nccreate(filename,param_props.file_name,'Dimensions',{'lon' xdim 'lat' ydim 'pres' zdim});
    ncwrite(filename,param_props.file_name,gobai_3d);

    % display information
    fprintf([alg_type ' Prediction (Month ' num2str(m) ', Week ' num2str(w) ')\n']);

end

%% for creating main netCDF file
function create_nc_file(TS,base_grid,xdim,ydim,zdim,gobai_alg_dir,...
    param_props)

% define file name
filename = [gobai_alg_dir 'gobai-' param_props.file_name '.nc'];

% create folder and file
if ~isfolder([gobai_alg_dir]); mkdir(gobai_alg_dir); end
if isfile(filename); delete(filename); end % delete file if it exists
% bgc parameter
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
else
    nccreate(filename,param_props.file_name,'Dimensions',{'lon',xdim,'lat',ydim,'depth',zdim,'time',Inf},...
        'DataType','single','FillValue',NaN);
end
ncwriteatt(filename,param_props.file_name,'units',param_props.units);
ncwriteatt(filename,param_props.file_name,'long_name',param_props.long_param_name);
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
ncwriteatt(filename,'time','units','days since 1950-0-0');
ncwriteatt(filename,'time','axis','T');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','_CoordinateAxisType','Time');

end

% train_ffnn
%
% DESCRIPTION:
% This function uses the combined dataset to train
% machine learning models in GMM clusters.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/17/2023

%% load combined data
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
load(['Data/processed_all_o2_data_' file_date float_file_ext '.mat'],...
     'all_data','file_date');

%% load data clusters
load(['Data/all_data_clusters_' num2str(num_clusters) '_' ...
    file_date float_file_ext '.mat'],'all_data_clusters');

%% remove float data for GLODAP only test
if glodap_only
    glodap_idx = all_data.platform < 10^6;
    vars = fieldnames(all_data);
    for v = 1:length(vars)
        all_data.(vars{v}) = all_data.(vars{v})(glodap_idx);
    end
end
clear glodap_only glodap_idx vars v

%% create directory and file names
ffnn_dir = ['Models/FFNN/FFNN_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio)];
ffnn_fnames = cell(num_clusters,1);
for c = 1:num_clusters
    ffnn_fnames(c) = ...
        {['FFNN_oxygen_C' num2str(c)]};
end
gobai_ffnn_dir = ...
    ['Data/GOBAI/FFNN_c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];

%% fit models using all data (NN)
% define model parameters
nodes1 = [5 10 15];
nodes2 = [15 10 5];

% fit models for each cluster
for c = 1:num_clusters

    % start timing fit
    tic

    % fit model for each cluster
    FFNN = ...
        fit_FFNN('oxygen',all_data,all_data_clusters.(['c' num2str(c)]),...
        true(size(all_data.platform)),variables,nodes1,nodes2,...
        train_ratio,val_ratio,test_ratio,thresh);

    % stop timing fit
    fprintf(['Train FFNN - Cluster #' num2str(c) ': ']);
    toc

    % save model for each cluster
    if ~isfolder([pwd '/' ffnn_dir]); mkdir(ffnn_dir);end
    save([ffnn_dir '/' ffnn_fnames{c}],'FFNN','-v7.3');

    % clean up
    clear FFNN

end

% clean up
clear all_data all_data_clusters

%% process time
start_year = 2004;
end_year = floor(snap_date/1e2);
end_month = mod(snap_date,1e2);
years = repelem(start_year:end_year,12)';
months = repmat(1:12,1,length(years)/12)';
years = years(1:end-(12-end_month));
months = months(1:end-(12-end_month));
clear start_year end_year end_month

%% load temperature and salinity data

% compute estimates for each month
parfor m = 1:length(years)

% import RG or RFROM
if strcmp(base_grid,'RG')

    % file names
    fname_temp = 'Data/RG_CLIM/RG_Climatology_Temp.nc';
    fname_sal = 'Data/RG_CLIM/RG_Climatology_Sal.nc';

    % load dimensions
    TS.lon = ncread(fname_temp,'Longitude');
    TS.lat = ncread(fname_temp,'Latitude');
    TS.pres = ncread(fname_temp,'Pressure');
    TS.time = ncread(fname_temp,'Time');

    % load variables
    TS.temp = ncread(fname_temp,'Temperature',[1 1 1 m],[Inf Inf Inf 1]);
    TS.sal = ncread(fname_sal,'Salinity',[1 1 1 m],[Inf Inf Inf 1]);

elseif strcmp(base_grid,'RFROM')

end

% define dimensions
xdim = length(TS.lon);
ydim = length(TS.lat);
zdim = length(TS.pres);

% replicate grid variables
TS.longitude = repmat(TS.lon,1,ydim,zdim);
TS.latitude = repmat(TS.lat',xdim,1,zdim);
TS.pressure = repmat(permute(TS.pres,[3 2 1]),xdim,ydim,1);

% convert to arrays
TS_index = ~isnan(TS.temp);
vars = fieldnames(TS);
for v = 1:length(vars)
    if ndims(TS.(vars{v})) == 3
        TS.([vars{v} '_array']) = TS.(vars{v})(TS_index);
        TS = rmfield(TS,vars{v});
    end
end

% process longitude
TS.longitude_array = convert_lon(TS.longitude_array);

% calculate absolute salinity, conservative temperature, potential density
TS.salinity_abs_array = gsw_SA_from_SP(TS.sal_array,TS.pressure_array,...
    TS.longitude_array,TS.latitude_array);
TS.temperature_cns_array = gsw_CT_from_t(TS.salinity_abs_array,...
    TS.temp_array,TS.pressure_array);
TS.sigma_array = gsw_sigma0(TS.salinity_abs_array,TS.temperature_cns_array);

% define variables for predictions
variables_TS = cell(size(variables));
for v = 1:length(variables)
    variables_TS{v} = [variables{v} '_array'];
end

%% apply models to gridded data (NN)

% pre-allocate
gobai_matrix = single(nan(length(TS.temp_array),num_clusters));
probs_matrix = single(nan(length(TS.temp_array),num_clusters));

% apply models for each cluster
for c = 1:num_clusters

    % load GMM cluster probabilities for this cluster and month, and convert to array
    GMM_probs = ...
        load(['Data/GMM_' num2str(num_clusters) '/c' num2str(c) '/m' num2str(m)],'GMM_monthly_probs');
    GMM_probs.monthly_probabilities = GMM_probs.GMM_monthly_probs;
    GMM_probs.probabilities_array = GMM_probs.monthly_probabilities(TS_index);
    GMM_probs = rmfield(GMM_probs,{'GMM_monthly_probs' 'monthly_probabilities'});
    probs_matrix(:,c) = GMM_probs.probabilities_array;

    % start timing predictions
    tic

    % load model for this cluster
    alg = load([ffnn_dir '/' ffnn_fnames{c}],'FFNN');

    % predict data for each cluster
    gobai_matrix(:,c) = ...
        run_FFNN(alg.FFNN,TS,GMM_probs.probabilities_array,...
        true(size(TS.temp_array)),variables_TS,thresh);

    % stop timing predictions
    fprintf(['Run FFNN - Cluster #' num2str(c) ': ']);
    toc

end

%% calculate weighted average over each cluster using probabilities

% weighted mean
gobai_array = ...
    double(sum(gobai_matrix.*probs_matrix,2,'omitnan')./...
    sum(probs_matrix,2,'omitnan'));

% back to 3D grid
gobai_3d = nan(length(TS.lon),length(TS.lat),length(TS.pres));
gobai_3d(TS_index) = gobai_array;

% save monthly output
if ~isfolder([pwd '/' gobai_ffnn_dir]); mkdir(gobai_ffnn_dir); end
parsave_v1([gobai_ffnn_dir 'm' num2str(m)],gobai_3d);

end

% end parallel session
delete(gcp('nocreate'));

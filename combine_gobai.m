% combine_gobai
%
% DESCRIPTION:
% This function combines output from gobai gridded
% fields obtained via different machine learning algorithms.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/22/2023

function combine_gobai(param,dir_base,base_grid,file_date,float_file_ext,...
    num_clusters,variables,train_ratio,val_ratio,test_ratio,thresh,...
    numWorkers_predict,years_to_predict)

% process parameter name
param1 = param_name(param);

%% create directory names
gobai_ffnn_dir = ... % FFNN
    [param1 '/Data/GOBAI/' base_grid '/FFNN/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/train' num2str(100*train_ratio) '_val' ...
    num2str(100*val_ratio) '_test' num2str(100*val_ratio) '/'];
gobai_rfr_dir = ... % RFR
    [param1 '/Data/GOBAI/' base_grid '/RFR/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numtrees) '_lf' num2str(minLeafSize) '/'];
gobai_gbm_dir = ... % GBM
    [param1 '/Data/GOBAI/' base_grid '/GBM/c' num2str(num_clusters) '_' file_date ...
    float_file_ext '/tr' num2str(numstumps) '/'];
gobai_dir = ... % final product
    [param1 '/Data/GOBAI/' base_grid '/' num2str(num_clusters) '_' file_date ...
    float_file_ext '/'];

%% process time
% start_year = 2004;
% end_year = floor(snap_date/1e2);
% end_month = mod(snap_date,1e2);
% years = repelem(start_year:end_year,12)';
% months = repmat(1:12,1,length(years)/12)';
% years = years(1:end-(12-end_month));
% months = months(1:end-(12-end_month));
% clear start_year end_year end_month

%% determine timesteps
if strcmp(base_grid,'RG')
    [~,timesteps] = load_RG_dim([pwd '/Data/RG_CLIM/']);
elseif strcmp(base_grid,'RFROM')
    [~,timesteps] = load_RFROM_dim([pwd '/Data/RFROM/']);
end

%% average over each time window
parfor m = 1:length(years)
    
    % load monthly output (FFNN)
    gobai_3d_ffnn = load([gobai_ffnn_dir 'm' num2str(m)],'x');
    % load monthly output (RFR)
    gobai_3d_rfr = load([gobai_rfr_dir 'm' num2str(m)],'x');
    % load monthly output (GBM)
    gobai_3d_gbm = load([gobai_gbm_dir 'm' num2str(m)],'x');

    % average monthly outputs
    gobai_3d_ens = mean(cat(4,gobai_3d_ffnn,gobai_3d_rfr,gobai_3d_gbm));

    % create folder and save monthly output
    if ~isfolder([pwd '/' gobai_dir]); mkdir(gobai_dir); end
    filename = [gobai_dir 'm' num2str(m) '_w' num2str(w) '.nc'];
    if isfile(filename); delete(filename); end % delete file if it exists
    % oxygen
    nccreate(filename,'o2','Dimensions',{'lon',xdim,'lat',ydim,'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'o2',gobai_3d);
    ncwriteatt(filename,'o2','units','umol/kg');
    ncwriteatt(filename,'o2','long_name','Dissolved Oxygen Amount Content');
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
    % pressure
    nccreate(filename,'pres','Dimensions',{'pres',zdim},...
        'DataType','single','FillValue',NaN);
    ncwrite(filename,'pres',TS.Pressure);
    ncwriteatt(filename,'pres','units','decibars');
    ncwriteatt(filename,'pres','axis','Z');
    ncwriteatt(filename,'pres','long_name','pressure');
    ncwriteatt(filename,'pres','_CoordinateAxisType','Pres');

end

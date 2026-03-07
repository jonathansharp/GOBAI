% process_gobai_netcdf
%
% DESCRIPTION:
% This function combines output from gobai gridded
% fields obtained via different machine learning algorithms.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 3/10/2025

function process_gobai_netcdf(alg_type,param_props,param_path,temp_path,sal_path,...
    base_grid,file_date,float_file_ext,num_clusters,variables,thresh,...
    numWorkers_predict,clust_vars,start_year,end_year,snap_date,varargin)

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

%% directory base
if strcmp(alg_type,'FFNN')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;train_ratio;val_ratio;test_ratio});
elseif strcmp(alg_type,'RFR')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
elseif strcmp(alg_type,'GBM')
    dir_base = create_dir_base(alg_type,{base_grid;num_clusters;file_date;...
        float_file_ext;numstumps;numbins});
end

%% establish paths, folders, and options
% file paths and names
fpath = [param_path dir_base];

% create annual file folder
mkdir([param_path dir_base '/gobai_annual']);
mkdir([param_path dir_base '/gobai_monthly']);

%% create monthly files
% loop through each year to create annual files
for y = start_year:end_year
    % pre-allocate annual monthly RFROMs
    
    % extract weekly salinity values and average to monthly
    for m = 1:12
        RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'ocean_salinity');
        RFROM.(['m' num2str(y)]) = ...
            cat(4,RFROM.(['m' num2str(y)]),mean(RFROM_temp,4,'omitnan'));
        RFROM_temp = ncread([fpath folder fname num2str(y) '_' ...
            sprintf('%02s',num2str(m)) '.nc'],'time');
        RFROM.(['m' num2str(y) '_time']) = ...
            [RFROM.(['m' num2str(y) '_time']);mean(RFROM_temp,'omitnan')];
        clear RFROM_temp
    end
    % copy schema to new annual file
    schema = ncinfo([fpath folder fname num2str(y) '_01.nc']);
    schema.Dimensions(3).Length = 12;
    schema.Variables(3).Size = 12;
    schema.Variables(3).Dimensions.Length = 12;
    schema.Variables(5).Size = [1440,720,58,12];
    schema.Variables(5).Dimensions(4).Length = 12;
    if isfile([fpath folder '_annual' fname num2str(y) '.nc'])
        delete([fpath folder '_annual' fname num2str(y) '.nc']);
    end
    ncwriteschema([fpath folder '_annual' fname num2str(y) '.nc'],schema);
    clear schema
    % add dimensional variables to monthly file
    vars = {'longitude' 'latitude' 'mean_pressure' 'mean_pressure_bnds'};
    for v = 1:length(vars)
        a=ncread([fpath folder fname num2str(y) '_01.nc'],vars{v});
        ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],vars{v},a);
    end
    % add averages to monthly file
    ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
        'time',RFROM.(['m' num2str(y) '_time']));
    ncwrite([fpath folder '_annual' fname num2str(y) '.nc'],...
        'ocean_salinity',RFROM.(['m' num2str(y)]));
    % clean up
    clear RFROM a m v y vars
end


% acquire_o2_snapshot_data
%
% DESCRIPTION:
% This function is used to import a monthly BGC Argo snapshot file, then
% extract all float data, filter the data by quality flag and data mode, do
% some additional QC and data removal, interpolate the data onto standard
% depth levels, and save and plot the processed float data file.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

function acquire_o2_snapshot_data(data_modes,float_file_ext,snap_date,snap_download)

%% Download BGC Argo snapshot
if snap_download == 1
    
    % set search patterns for parsing of web page
    pattern = 'BGC Sprof data files (';
    pat_month = 'BGC Sprof data files \(20\d{2}-\d{2}-\d{2} snapshot\)';
    mth_idx = length(pattern) + 1; % first character of YYYY-MM
    
    % download the web page that contains the links to all snapshots
    options = weboptions; options.Timeout = 10; % increase timeout options
    page = webread('https://www.seanoe.org/data/00311/42182/');
    idx = strfind(page, pattern);
    if isempty(idx)
        fprintf('The snapshot web page could not be read properly.\n')
        return
    end
    clear pattern options 
    
    % determine which snapshots are available
    snap_month = [];
    snap_size = [];
    snap_url = {};
    count = 0;
    for i = 1:length(idx)
        this_link = page(idx(i):idx(i)+400);
        match_month = regexp(this_link, pat_month, 'match', 'once');
        pat_url = '"fileUrl":"http[\w/:.]+"';
        match_url = regexp(this_link, pat_url, 'match', 'once');
        pat_size = '"size":\d+';
        match_size = regexp(this_link, pat_size, 'match', 'once');
        % note that some snapshots are only available on demand, for those,
        % there is no fileUrl entry
        if ~isempty(match_month) && ~isempty(match_url)
            count = count + 1;
            year = str2double(match_month(mth_idx:mth_idx+3));
            month = str2double(match_month(mth_idx+5:mth_idx+6));
            snap_month(count) = 100 * year + month;
            snap_url{count} = match_url(12:end-1);
            if ~isempty(match_size)
                snap_size(count) = uint64(str2double(match_size(8:end)));
            else
                snap_size(count) = -1;
            end
        end
    end
    clear pat_month mth_idx page idx snap_size this_link match_month
    clear pat_url match_url pat_size match_size year month i
    
    % return error if no snapshots are found
    if ~count
        fprintf('No matching snapshots were found.\n')
        return
    end
    
    % use the most recent snapshot, need to sort them by snap_date
    [~,isort] = sort(snap_month,'descend'); % most recent will be first
    isnap = isort(1);
    snap_date = snap_month(isnap); % assignment makes it easier going forward
    clear snap_month count isort
    
    % this is the most commonly used format for file names:
    filename = regexp(snap_url{isnap}, '\d+\.tar.gz', 'match', 'once');
    if isempty(filename) % alternate file name format
        filename = regexp(snap_url{isnap}, '\d+\.tgz', 'match', 'once');
    end
    
    % Download snapshot from online DOI if needed
    if ~exist('BGC_Argo_Snapshots','dir'); mkdir('BGC_Argo_Snapshots'); end
    if exist(['BGC_Argo_Snapshots/' num2str(snap_date) '-BgcArgoSprof/'],'dir') ~= 7
        url = 'https://www.seanoe.org/data/00311/42182/data/';
        websave(filename,[url filename]); % doenload zipped file
        gunzip(filename,'BGC_Argo_Snapshots'); % unzip file
        % define tarball name
        tarname = regexp(filename, '\d+\.tar', 'match', 'once');
        if isempty(filename) % alternate file name format
            filename = regexp(snap_url{isnap}, '\d+\.tgz', 'match', 'once');
        end
        % untar tarball
        untar(['BGC_Argo_Snapshots/' tarname],'BGC_Argo_Snapshots');
        % clean up
        delete(filename);
        delete(['BGC_Argo_Snapshots/' tarname]);
    end
    clear url snap_url isnap filename tarname
else
    snap_date = snap_date;
end

%% Only do all this if downloaded float matlab file does not exist
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
% if exist(['O2/Data/processed_float_o2_data_' file_date float_file_ext '.mat'],'file') ~= 2

%% Define file structure
snapshot_path = ['BGC_Argo_Snapshots/' num2str(snap_date) '-BgcArgoSprof/dac']; % define file path
folderinfo = dir(snapshot_path); % extract info about folders and files
foldernames = {folderinfo.name}'; % define folder and file names
idx_folders = find(~contains(foldernames,'.')); % index to folders only
clear today_date mnth folderinfo

%% pre-allocate float data structure
float_data.OXY = [];
float_data.OXY_LAT = [];
float_data.OXY_LON = [];
float_data.OXY_PRES = [];
float_data.OXY_TIME = [];
float_data.OXY_TEMP = [];
float_data.OXY_SAL = [];
float_data.OXY_FLOAT = [];
float_data.OXY_PROF_ID = [];

%% set up figure
figure(1); hold on;
set(gcf,'visible','off','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch','k');
m_grid('linestyle','-','ytick',-90:30:90);

%% Define interpolation parameters
vars = {'PSAL' 'TEMP' 'DOXY'}; % define variables to interpolate
zi = ([2.5 10:10:170 182.5 ... % construct depth axis on which to interpolate
    200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';

%% Download and process all float data
for n = 1:length(idx_folders) % for each DAC

    %% Define file info
    % extract float file info
    fileinfo = dir([snapshot_path '/' foldernames{idx_folders(n)} '/*.nc']);
    filenames = {fileinfo.name}'; % define Sprof file names
    clear fileinfo

    %% Download and process floats in each DAC
    for f = 1:length(filenames) % for each float

        %% Download float data
        floatnum = extractBefore(filenames{f},'_'); % define float number
        float.(['F' floatnum]) = ... % download float data
            netcdfreader([snapshot_path '/' foldernames{idx_folders(n)} '/' filenames{f}]);
        [n_prof, n_param, n_levels] = ... % get dimensions
            get_dims([snapshot_path '/' foldernames{idx_folders(n)} '/' filenames{f}]);

        %% Parse parameter names
        temp = cell(n_param, 1);
        % extract parameter names as coherent strings
        for m = 1:n_param
            temp{m} = strrep(float.(['F' floatnum]).PARAMETER(:,m)', ' ', '');
        end
        float.(['F' floatnum]).PARAMETER = temp;
        clear temp m
        
        %% Continue if float has an oxygen sensor
        if any(strcmp(float.(['F' floatnum]).PARAMETER,'DOXY'))

            %% Parse parameter data modes
            for m = 1:n_param
                if ~isempty(float.(['F' floatnum]).PARAMETER{m})
                    float.(['F' floatnum]).([float.(['F' floatnum]).PARAMETER{m} '_DATA_MODE']) = ...
                        repmat(float.(['F' floatnum]).('PARAMETER_DATA_MODE')(m,:),n_levels,1);
                end
            end
            clear m
    
            %% Convert quality flag strings to numeric format
            all_vars = fieldnames(float.(['F' floatnum]));
            for v = 1:length(all_vars)
                if endsWith(all_vars{v},'_QC') && ... % Check for QC identifier
                        ~startsWith(all_vars{v},'PROFILE') % But not a profile QC
                    if isequal(size(float.(['F' floatnum]).(all_vars{v})), [n_prof 1])
                        float.(['F' floatnum]).(all_vars{v}) = ...
                            repmat(float.(['F' floatnum]).(all_vars{v})',n_levels,1);
                    end
                    float.(['F' floatnum]).(all_vars{v}) = ...
                        float.(['F' floatnum]).(all_vars{v})(:);
                    float.(['F' floatnum]).(all_vars{v}) = ...
                        strrep(float.(['F' floatnum]).(all_vars{v})', ' ', '0')';
                    float.(['F' floatnum]).(all_vars{v}) = ...
                        str2num(float.(['F' floatnum]).(all_vars{v}));
                    float.(['F' floatnum]).(all_vars{v}) = ...
                        reshape(float.(['F' floatnum]).(all_vars{v}), n_levels, n_prof);
                end
            end
            clear all_vars v

            %% add unprocessed oxygen data to plot
            lon_temp = convert_lon(convert_lon(float.(['F' floatnum]).LONGITUDE));
            m_scatter(lon_temp,float.(['F' floatnum]).LATITUDE,'.k');
            clear lon_temp

            %% Loop trough floats, interpolate profiles, and log data
            for k = 1:numel(vars) % for each variable
                % index based on data mode
                mode_idx = false(n_levels,n_prof);
                for m = 1:length(data_modes)
                    mode_idx = mode_idx | ...
                        float.(['F' floatnum]).([vars{k} '_DATA_MODE']) == data_modes{m};
                end
                % index based on data mode and flags
                index = mode_idx & ... % use data mode index that was created
                        (float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 1 | ...
                        float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 2 | ...
                        float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 8) & ... % with good variable flags,
                        (float.(['F' floatnum]).PRES_ADJUSTED_QC == 1 | ...
                        float.(['F' floatnum]).PRES_ADJUSTED_QC == 2); % and good pressure flags.
                
                % pre-allocate interpolated profiles
                float.(['F' floatnum]).([vars{k} '_ADJUSTEDi']) = ...
                        nan(length(zi),length(float.(['F' floatnum]).CYCLE_NUMBER));
                % pre-allocate nan profile index
                nan_idx = false(length(float.(['F' floatnum]).CYCLE_NUMBER),1);
                
                % loop through profiles and interpolate
                for p = 1:length(float.(['F' floatnum]).CYCLE_NUMBER) % for each profile
                    if sum(index(:,p)) > 10 % if more than ten data points are available
                        % interpolate to edges (without extrapolation)
                        temp_pres = float.(['F' floatnum]).PRES_ADJUSTED(index(:,p),p);
                        temp_var = float.(['F' floatnum]).([vars{k} '_ADJUSTED'])(index(:,p),p);
                        [~,unique_idx_pres] = unique(temp_pres);
                        temp_var_i = interp1(temp_pres(unique_idx_pres),...
                            temp_var(unique_idx_pres),zi,'linear');
        
                        % log interpolated profile
                        float.(['F' floatnum]).([vars{k} '_ADJUSTEDi'])(:,p) = temp_var_i;
        
                    else
        
                        % log profile as not available if no matching data
                        nan_idx(p) = true;
                        
                    end
                end
            end
            clear index p k temp_pres temp_var unique_idx_pres temp_var_i
    
            %% calculate time from JULD
            float.(['F' floatnum]).TIME = datenum(1950,1,1+float.(['F' floatnum]).JULD);
    
            %% drop empty interpolated profiles
            float.(['F' floatnum]).DOXY_ADJUSTEDi(:,nan_idx) = [];
            float.(['F' floatnum]).PSAL_ADJUSTEDi(:,nan_idx) = [];
            float.(['F' floatnum]).TEMP_ADJUSTEDi(:,nan_idx) = [];
    
            %% add other variables to replicate interpolated variable size
            float.(['F' floatnum]).PRES_ADJUSTEDi = repmat(zi,1,sum(~nan_idx));
            float.(['F' floatnum]).LATITUDEi = repmat(float.(['F' floatnum]).LATITUDE(~nan_idx)',length(zi),1);
            float.(['F' floatnum]).LONGITUDEi = repmat(float.(['F' floatnum]).LONGITUDE(~nan_idx)',length(zi),1);
            float.(['F' floatnum]).TIMEi = repmat(float.(['F' floatnum]).TIME(~nan_idx)',length(zi),1);
            float.(['F' floatnum]).CYCLE_NUMBERi = repmat(float.(['F' floatnum]).CYCLE_NUMBER(~nan_idx)',length(zi),1);
            float.(['F' floatnum]).FLOAT_NUMBERi = repmat(str2double(floatnum),length(zi),sum(~nan_idx));
    
            %% add interpolated data to float data structure
            float_data.OXY = [float_data.OXY;float.(['F' floatnum]).DOXY_ADJUSTEDi(:)];
            float_data.OXY_LAT = [float_data.OXY_LAT;float.(['F' floatnum]).LATITUDEi(:)];
            float_data.OXY_LON = [float_data.OXY_LON;float.(['F' floatnum]).LONGITUDEi(:)];
            float_data.OXY_PRES = [float_data.OXY_PRES;float.(['F' floatnum]).PRES_ADJUSTEDi(:)];
            float_data.OXY_TIME = [float_data.OXY_TIME;float.(['F' floatnum]).TIMEi(:)];
            float_data.OXY_TEMP = [float_data.OXY_TEMP;float.(['F' floatnum]).TEMP_ADJUSTEDi(:)];
            float_data.OXY_SAL = [float_data.OXY_SAL;float.(['F' floatnum]).PSAL_ADJUSTEDi(:)];
            float_data.OXY_FLOAT = [float_data.OXY_FLOAT;float.(['F' floatnum]).FLOAT_NUMBERi(:)];
            float_data.OXY_PROF_ID = [float_data.OXY_PROF_ID;float.(['F' floatnum]).FLOAT_NUMBERi(:).*1000 + ...
                float.(['F' floatnum]).CYCLE_NUMBERi(:)];

        end

        %% clean up
        clear float floatnum nan_idx n_prof n_param n_levels

    end

    %% clean up
    clear f filenames

end

%% clean up
clear foldernames snapshot_path idx_folders n vars zi

%% remove nan data points
idx_oxy = isnan(float_data.OXY);
idx_sal = isnan(float_data.OXY_SAL);
idx_temp = isnan(float_data.OXY_TEMP);
idx_any = idx_oxy | idx_sal | idx_temp;
vars = fieldnames(float_data);
for v = 1:length(vars)
    float_data.(vars{v})(idx_any) = [];
end
clear idx_oxy idx_sal idx_temp idx_any v vars

%% plot processed oxygen data
figure(1); hold on;
lon_temp = convert_lon(convert_lon(float_data.OXY_LON));
m_scatter(lon_temp,float_data.OXY_LAT,'.g');
if ~exist('Figures','dir'); mkdir('Figures'); end
exportgraphics(gcf,['O2/Figures/Data/processed_float_'  file_date float_file_ext '.png']);
close
clear lon_temp

%% convert times and dates
date_temp = datevec(float_data.OXY_TIME);
date_temp0 = date_temp;
date_temp0(:,2:3) = 1; % Jan. 1 of each year
float_data.OXY_YEAR = date_temp(:,1);
float_data.OXY_DAY = datenum(date_temp) - datenum(date_temp0) + 1;
clear date_temp date_temp0

%% Calculate absolute salinity, conservative temperature, potential density, and spice
float_data.OXY_ABSSAL = gsw_SA_from_SP(float_data.OXY_SAL,float_data.OXY_PRES,float_data.OXY_LON,float_data.OXY_LAT);
float_data.OXY_CNSTEMP = gsw_CT_from_t(float_data.OXY_ABSSAL,float_data.OXY_TEMP,float_data.OXY_PRES);
float_data.OXY_SIGMA = gsw_sigma0(float_data.OXY_ABSSAL,float_data.OXY_CNSTEMP);
float_data.OXY_SPICE = gsw_spiciness0(float_data.OXY_ABSSAL,float_data.OXY_CNSTEMP);

%% display the number of matching floats and profiles
disp(['# of matching Argo profiles (OXY): ' num2str(length(unique(float_data.OXY_PROF_ID)))]);
disp(['# of matching Argo floats (OXY): ' num2str(length(unique(float_data.OXY_FLOAT)))]);

%% save processed float data
if ~exist([pwd '/O2/Data'],'dir'); mkdir('O2/Data'); end
save(['O2/Data/processed_float_o2_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date','-v7.3');

%% clean up
clear float_data
close all

% end

end
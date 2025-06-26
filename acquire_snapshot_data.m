% acquire_snapshot_data
%
% DESCRIPTION:
% This function is used to import a monthly BGC Argo snapshot file, then
% extract all float data, filter the data by quality flag and data mode, do
% some additional QC and data removal, interpolate the data onto standard
% depth levels, and save and plot the processed float data file.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 4/11/2025

function acquire_snapshot_data(param_props,data_modes,float_file_ext,snap_date,snap_download)

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

    % use the year/month snapshot input
    isnap = find(snap_month == snap_date);
    clear snap_month count

%     % use the most recent snapshot, need to sort them by snap_date
%     [~,isort] = sort(snap_month,'descend'); % most recent will be first
%     isnap = isort(1);
%     snap_date = snap_month(isnap); % assignment makes it easier going forward
%     clear snap_month count isort
    
    % this is the most commonly used format for file names:
    filename = regexp(snap_url{isnap}, '\d+\.tar.gz', 'match', 'once');
    if isempty(filename) % alternate file name format
        filename = regexp(snap_url{isnap}, '\d+\.tgz', 'match', 'once');
    end
    
    % Download snapshot from online DOI if needed
    if ~exist('BGC_Argo_Snapshots','dir'); mkdir('BGC_Argo_Snapshots'); end
    if exist(['BGC_Argo_Snapshots/' num2str(snap_date) '-BgcArgoSprof/'],'dir') ~= 7
        url = 'https://www.seanoe.org/data/00311/42182/data/';
        websave(filename,[url filename]); % download zipped file
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
        disp('Float snapshot acquired.')
    else
        disp('Float snapshot already downloaded.')
    end
    clear url snap_url isnap filename tarname
    
else
    disp('Float snapshot download not requested.')
end

%% Only do all this if downloaded float matlab file does not exist
file_date = datestr(datenum(floor(snap_date/1e2),mod(snap_date,1e2),1),'mmm-yyyy');
if exist([param_props.dir_name '/Data/processed_float_' param_props.file_name '_data_' file_date float_file_ext '.mat'],'file') ~= 2

%% Define file structure
snapshot_path = ['BGC_Argo_Snapshots/' num2str(snap_date) '-BgcArgoSprof/dac']; % define file path
folderinfo = dir(snapshot_path); % extract info about folders and files
foldernames = {folderinfo.name}'; % define folder and file names
idx_folders = find(~contains(foldernames,'.')); % index to folders only
clear today_date mnth folderinfo

%% pre-allocate float data structure
if strcmp(param_props.argo_name,'PH_IN_SITU_TOTAL')
    float_data.OXY = [];
elseif strcmp(param_props.argo_name,'NITRATE')
    float_data.OXY = [];
end
float_data.(param_props.temp_name) = [];
float_data.LAT = [];
float_data.LON = [];
float_data.PRES = [];
float_data.TIME = [];
float_data.TEMP = [];
float_data.SAL = [];
float_data.FLOAT = [];
float_data.PROF_ID = [];

%% set up figure
figure(1); hold on;
set(gcf,'visible','off','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);

%% Define interpolation parameters
if strcmp(param_props.argo_name,'PH_IN_SITU_TOTAL')
    vars = {'PSAL' 'TEMP' 'DOXY' param_props.argo_name}; % define variables to interpolate
elseif strcmp(param_props.argo_name,'NITRATE')
    vars = {'PSAL' 'TEMP' 'DOXY' param_props.argo_name}; % define variables to interpolate
else
    vars = {'PSAL' 'TEMP' param_props.argo_name}; % define variables to interpolate
end
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
    counter = 1;
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
        
        %% Continue if float has the relevant sensors
        cmp_idx = true;
        for cmp_n = 1:length(vars)
            % if any needed parameter is not present, change index to false
            if any(strcmp(float.(['F' floatnum]).PARAMETER,vars(cmp_n)))
            else; cmp_idx = false; end
        end

        %% Extract float data
        if cmp_idx
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

            %% add unprocessed data to plot
            lon_temp = convert_lon(convert_lon(float.(['F' floatnum]).LONGITUDE));
            lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
            m_scatter(lon_temp,float.(['F' floatnum]).LATITUDE,'.k');
            clear lon_temp

            %% Loop through floats, interpolate profiles, and log data
            for k = 1:numel(vars) % for each variable
                % index based on data mode
                mode_idx = false(n_levels,n_prof);
                for m = 1:length(data_modes)
                    if ~isfield(float.(['F' floatnum]),([vars{k} '_DATA_MODE']))
                        mode_idx = mode_idx;
                    else
                    mode_idx = mode_idx | ...
                        float.(['F' floatnum]).([vars{k} '_DATA_MODE']) == data_modes{m};
                    end
                end
                % index based on data mode and flags and depth < 1950
                index = mode_idx & ... % use data mode index that was created
                        (float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 1 | ...
                        float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 2 | ...
                        float.(['F' floatnum]).([vars{k} '_ADJUSTED_QC']) == 8) & ... % with good variable flags,
                        (float.(['F' floatnum]).PRES_ADJUSTED_QC == 1 | ...
                        float.(['F' floatnum]).PRES_ADJUSTED_QC == 2) & ... % and good pressure flags
                        float.(['F' floatnum]).PRES_ADJUSTED <= 1950;
                
                % pre-allocate interpolated profiles
                float.(['F' floatnum]).([vars{k} '_ADJUSTEDi']) = ...
                        nan(length(zi),length(float.(['F' floatnum]).CYCLE_NUMBER));
                % pre-allocate nan profile index
                nan_idx = false(length(float.(['F' floatnum]).CYCLE_NUMBER),1);
                
                % loop through profiles and interpolate
                for p = 1:length(float.(['F' floatnum]).CYCLE_NUMBER) % for each profile
                    if sum(index(:,p)) > 10 % if more than ten data points are available
                        % extract temporary pressure axis and parameter
                        temp_pres = float.(['F' floatnum]).PRES_ADJUSTED(index(:,p),p);
                        temp_var = float.(['F' floatnum]).([vars{k} '_ADJUSTED'])(index(:,p),p);

                        % interpolate to edges (with extrapolation)
                        [~,unique_idx_pres] = unique(temp_pres);
                        try
                        temp_var_i = interp1(temp_pres(unique_idx_pres),...
                            temp_var(unique_idx_pres),zi,'linear','extrap');
                        catch
                        keyboard
                        end

                        % remove interpolated data beneath deepest measurement
                        % nah this doesn't seem necessary
                        % temp_var_i(zi>max(temp_pres(unique_idx_pres))) = NaN;

                        % remove interpolated data more than 100m from a measurement
                        index_from = false(size(zi));
                        for z = 1:length(zi)
                            index_from(z) = any(abs(temp_pres-zi(z)) < 100);
                        end
                        temp_var_i(~index_from) = NaN;

                        % % not sure why this is here...
                        % if str2double(floatnum)*1000 + float.(['F' floatnum]).CYCLE_NUMBER(p) == 5901464005
                        %     keyboard
                        % end

                        % if there is a greater than 0.5 % change per meter
                        % at the bottom of the profile or 1% change per
                        % meter at the top, remove extrapolated values
                        pres_axis = temp_pres(unique_idx_pres);
                        var_axis = temp_var(unique_idx_pres);
                        % bottom
                        bot_gradient = abs((100*((var_axis(end)-var_axis(end-1))./...
                            var_axis(end)))./(pres_axis(end)-pres_axis(end-1)));
                        if bot_gradient > 0.5
                            temp_var_i(zi>max(pres_axis)) = NaN;
                        end
                        % top
                        top_gradient = abs((100*((var_axis(2)-var_axis(1))./...
                            var_axis(2)))./(pres_axis(2)-pres_axis(1)));
                        if top_gradient > 5
                            temp_var_i(zi<min(pres_axis)) = NaN;
                        end

                        % log interpolated profile
                        float.(['F' floatnum]).([vars{k} '_ADJUSTEDi'])(:,p) = temp_var_i;

                    else
        
                        % log profile as not available if no matching data
                        nan_idx(p) = true;
                        
                    end
                end

                % save a random profile every 100 floats
                if k == 3 && mod(counter,100) == 0
                    prof_to_plot = randi(n_prof);
                    figure('visible','off'); hold on; set(gca,'YDir','reverse');
                    plot(float.(['F' floatnum]).([vars{k} '_ADJUSTEDi'])(:,prof_to_plot),...
                        zi,'linewidth',2);
                    info = gcf; info.Position(4) = info.Position(4)*2;
                    scatter(float.(['F' floatnum]).([vars{k} '_ADJUSTED'])(:,prof_to_plot),...
                        float.(['F' floatnum]).PRES_ADJUSTED(:,prof_to_plot),'filled');
                    legend({'Interpolation' 'Measurements'},'Location',...
                        'northoutside','NumColumns',2);
                    if ~exist([param_props.dir_name '/Figures/Data/Profiles'],'dir')
                        mkdir([param_props.dir_name '/Figures/Data/Profiles']); end
                    export_fig(gcf,[param_props.dir_name '/Figures/Data/Profiles/float_' ...
                        num2str(floatnum) '_prof' num2str(prof_to_plot) '.png'],'-transparent');
                    close
                end
 
            end
            clear index p k temp_pres temp_var unique_idx_pres temp_var_i
    
            %% calculate time from JULD
            float.(['F' floatnum]).TIME = datenum(1950,1,1+float.(['F' floatnum]).JULD);
    
            %% drop empty interpolated profiles
            if strcmp(param_props.argo_name,'PH_IN_SITU_TOTAL')
                float.(['F' floatnum]).DOXY_ADJUSTEDi(:,nan_idx) = [];
            elseif strcmp(param_props.argo_name,'NITRATE')
                float.(['F' floatnum]).DOXY_ADJUSTEDi(:,nan_idx) = [];
            end
            float.(['F' floatnum]).([param_props.argo_name '_ADJUSTEDi'])(:,nan_idx) = [];
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
            if strcmp(param_props.argo_name,'PH_IN_SITU_TOTAL')
                float_data.OXY = [float_data.OXY;float.(['F' floatnum]).DOXY_ADJUSTEDi(:)];
            elseif strcmp(param_props.argo_name,'NITRATE')
                float_data.OXY = [float_data.OXY;float.(['F' floatnum]).DOXY_ADJUSTEDi(:)];
            end
            float_data.(param_props.temp_name) = [float_data.(param_props.temp_name);float.(['F' floatnum]).([param_props.argo_name '_ADJUSTEDi'])(:)];
            float_data.LAT = [float_data.LAT;float.(['F' floatnum]).LATITUDEi(:)];
            float_data.LON = [float_data.LON;float.(['F' floatnum]).LONGITUDEi(:)];
            float_data.PRES = [float_data.PRES;float.(['F' floatnum]).PRES_ADJUSTEDi(:)];
            float_data.TIME = [float_data.TIME;float.(['F' floatnum]).TIMEi(:)];
            float_data.TEMP = [float_data.TEMP;float.(['F' floatnum]).TEMP_ADJUSTEDi(:)];
            float_data.SAL = [float_data.SAL;float.(['F' floatnum]).PSAL_ADJUSTEDi(:)];
            float_data.FLOAT = [float_data.FLOAT;float.(['F' floatnum]).FLOAT_NUMBERi(:)];
            float_data.PROF_ID = [float_data.PROF_ID;float.(['F' floatnum]).FLOAT_NUMBERi(:).*1000 + ...
                float.(['F' floatnum]).CYCLE_NUMBERi(:)];

        end

    % clean up
    float = rmfield(float,(['F' floatnum]));

    % increase counter
    counter = counter + 1;

    end

    %% clean up
    clear f filenames

end

%% clean up
clear foldernames snapshot_path idx_folders n vars zi

%% remove nan data points
if strcmp(param_props.argo_name,'NITRATE')
    idx = isnan(float_data.OXY) & isnan(float_data.(param_props.temp_name));
else
    idx = isnan(float_data.(param_props.temp_name));
end
idx_sal = isnan(float_data.SAL);
idx_temp = isnan(float_data.TEMP);
idx_any = idx | idx_sal | idx_temp;
vars = fieldnames(float_data);
for v = 1:length(vars)
    float_data.(vars{v})(idx_any) = [];
end
clear idx idx_sal idx_temp idx_any v vars

%% plot processed data
figure(1); hold on;
lon_temp = convert_lon(convert_lon(float_data.LON));
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,float_data.LAT,'.g');
if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir')
    mkdir([param_props.dir_name '/Figures/Data']); end
export_fig(gcf,[param_props.dir_name '/Figures/Data/processed_float_'  ...
    file_date float_file_ext '.png'],'-transparent');
close
clear lon_temp

%% convert times and dates
date_temp = datevec(float_data.TIME);
date_temp0 = date_temp;
date_temp0(:,2:3) = 1; % Jan. 1 of each year
float_data.YEAR = date_temp(:,1);
float_data.DAY = datenum(date_temp) - datenum(date_temp0) + 1;
clear date_temp date_temp0

%% Calculate absolute salinity, conservative temperature, potential density, and spice
float_data.ABSSAL = gsw_SA_from_SP(float_data.SAL,float_data.PRES,float_data.LON,float_data.LAT);
float_data.CNSTEMP = gsw_CT_from_t(float_data.ABSSAL,float_data.TEMP,float_data.PRES);
float_data.SIGMA = gsw_sigma0(float_data.ABSSAL,float_data.CNSTEMP);
float_data.SPICE = gsw_spiciness0(float_data.ABSSAL,float_data.CNSTEMP);

%% display the number of matching floats and profiles
disp(['# of matching Argo profiles (' param_props.argo_name '): ' num2str(length(unique(float_data.PROF_ID)))]);
disp(['# of matching Argo floats (' param_props.argo_name '): ' num2str(length(unique(float_data.FLOAT)))]);

%% save processed float data
if ~exist([pwd '/' param_props.dir_name '/Data'],'dir')
    mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/processed_float_' param_props.file_name '_data_' file_date float_file_ext '.mat'],...
    'float_data','file_date','-v7.3');

%% clean up
clear float_data
close all

% display information
disp('Float data processed and saved.')

else

% display information
disp('Float data already processed.')

end

end
function wod = acquire_wod_data(param_props,glodap_year,start_year,osd_opt,ctd_opt)

%% process parameter name
if strcmp(param_props.file_name,'o2')
    var_name = 'Oxygen';
elseif strcmp(param_props.file_name,'no3')
    var_name = 'Nitrate';
end

%% load glodap data
load([param_props.dir_name '/Data/processed_glodap_' param_props.file_name '_data_' num2str(glodap_year) '.mat'],...
    'glodap_data');

%% Only do all this if downloaded wod matlab file does not exist
todays_date = datevec(date);
end_year = todays_date(1);
year = num2str(end_year);

% if exist([param_props.dir_name '/Data/processed_wod_osd_' param_props.file_name '_data_' year '.mat'],'file') ~= 2

%% load WOD OSD profile data
folder = 'Data/WOD_Profiles_data';
if osd_opt == 1 && ctd_opt == 1
    types = {'osd' 'ctd'};
elseif osd_opt == 1 && ctd_opt == 0
    types = {'osd'};
elseif osd_opt == 0 && ctd_opt == 1
    types = {'ctd'};
else
    types = {};
end

%% define variables
dims = {'wod_unique_cast' 'time' 'date' 'lat' 'lon'};
vars = {'Temperature' 'Salinity' 'z' var_name};
vars_both = [dims vars];

%% pre-allocate osd data structure
for v = 1:length(vars_both); wod.(vars_both{v}) = []; end
for v = 1:length(vars); wod.([vars{v} '_flag']) = []; end
wod_data.(param_props.temp_name) = [];
wod_data.LAT = [];
wod_data.LON = [];
wod_data.PRES = [];
wod_data.TIME = [];
wod_data.YEAR = [];
wod_data.DAY = [];
wod_data.TEMP = [];
wod_data.SAL = [];
wod_data.CRU = [];
wod_data.ID = [];

%% set up figure
figure(1); hold on;
set(gcf,'visible','off','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);

%% load variables
for y = start_year:end_year
    for x = 1:length(types)
        file = [folder '/wod_' types{x} '_' num2str(y) '.nc'];
        if ~exist(file,'file')
            try
                websave(file,['https://www.ncei.noaa.gov/data/oceans/ncei/wod/' ...
                    num2str(y) '/wod_' types{x} '_' num2str(y) '.nc']);
            catch
                disp(['No wod ' types{x} ' file for ' num2str(y)])
            end
        end
        if exist(file,'file')
            schema = ncinfo(file);
            pdim = schema.Dimensions(1).Length;
            zdim = schema.Dimensions(2).Length;
            % dimensions
            for d = 1:length(dims)
                wod_temp.(dims{d}) = ncread(file,dims{d});
            end
            
            % download variabless
            for v = 1:length(vars)
                wod_temp.(vars{v}) = ncread(file,vars{v});
                wod_temp.([vars{v} '_flag']) = ncread(file,[vars{v} '_WODflag']);
                wod_temp.([vars{v} '_rows']) = ncread(file,[vars{v} '_row_size']);
                wod_temp.([vars{v} '_rows'])(isnan(wod_temp.([vars{v} '_rows']))) = 0;
                wod_temp.([vars{v} '_casts']) = ...
                    wod_temp.wod_unique_cast(wod_temp.([vars{v} '_rows']) > 0);
            end

            % index to common casts
            if strcmp(param_props.file_name,'o2')
                idx = ismember(wod_temp.wod_unique_cast,wod_temp.Oxygen_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.Temperature_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.Salinity_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.z_casts);
            elseif strcmp(param_props.file_name,'no3')
                idx = ismember(wod_temp.wod_unique_cast,wod_temp.Nitrate_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.Oxygen_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.Temperature_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.Salinity_casts) & ...
                      ismember(wod_temp.wod_unique_cast,wod_temp.z_casts);
            end

            % pre-allocate data for common casts
            for v = 1:length(vars)
                wod_temp.([vars{v} '_idx']) = 1;
                wod_temp.([vars{v} '_common']) = [];
                wod_temp.([vars{v} '_flag_common']) = [];
            end
            for d = 1:length(dims)
                wod_temp.([dims{d} '_common']) = [];
            end

            % extract data from common casts
            for c = 1:length(wod_temp.wod_unique_cast)
                if idx(c) % check for common cast
                    for v = 1:length(vars)
                        % extract variables and flags for common casts
                        wod_temp.([vars{v} '_common']) = [wod_temp.([vars{v} '_common']);...
                            wod_temp.([vars{v}])(wod_temp.([vars{v} '_idx']):wod_temp.([vars{v} '_idx'])+wod_temp.([vars{v} '_rows'])(c)-1,:)];
                        wod_temp.([vars{v} '_flag_common']) = [wod_temp.([vars{v} '_flag_common']);...
                            wod_temp.([vars{v} '_flag'])(wod_temp.([vars{v} '_idx']):wod_temp.([vars{v} '_idx'])+wod_temp.([vars{v} '_rows'])(c)-1,:)];
                    end
                    % add dimensions for unique casts
                    for d = 1:length(dims)
                        wod_temp.([dims{d} '_common']) = [wod_temp.([dims{d} '_common']);...
                            repmat(wod_temp.([dims{d}])(c),wod_temp.([vars{v} '_rows'])(c),1)];
                    end
                end
                % increase index for each variable
                for v = 1:length(vars)
                    wod_temp.([vars{v} '_idx']) = wod_temp.([vars{v} '_idx'])+wod_temp.([vars{v} '_rows'])(c);
                end
            end

            % index to where all data points have a value and good flag
            idx = ~isnan(wod_temp.Temperature_common) & wod_temp.Temperature_flag_common == 0 & ...
                ~isnan(wod_temp.Salinity_common) & wod_temp.Salinity_flag_common == 0 & ...
                ~isnan(wod_temp.z_common) & wod_temp.z_flag_common == 0 & ...
                ~isnan(wod_temp.Oxygen_common) & wod_temp.Oxygen_flag_common == 0;
            for v = 1:length(vars)
                wod.(vars{v}) = [wod.(vars{v});...
                    wod_temp.([vars{v} '_common'])(idx)];
                wod.([vars{v} '_flag']) = [wod.([vars{v} '_flag']);...
                    wod_temp.([vars{v} '_flag_common'])(idx)];
            end
            for d = 1:length(dims)
                wod.(dims{d}) = [wod.(dims{d});...
                    wod_temp.([dims{d} '_common'])(idx)];
            end

        end
        disp([num2str(y) ' processed - ' types{x} '.'])
    end
end

keyboard

%% adjust time
wod.time = datenum(1770,1,1+wod.time);

%% calculate pressure from depth
wod.pres = gsw_p_from_z(-wod.z,wod.lat);

%% plot unprocessed profile locations
figure(1); hold on;
set(gcf,'visible','on','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
lon_temp = convert_lon(convert_lon(wod.lon,'0-360'));
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,wod.lat,'.k');
hold off;
clear lon_temp

%% construct depth axis on which to interpolate
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
% define variables to interpolate
vars = {var_name 'Salinity' 'Temperature'};
varsi = {param_props.temp_name 'SAL' 'TEMP'};
% define station ids
stations = unique(wod.profile);

%% process osd data
for f = 1:length(stations) % for each unique station id

    %% index according to station id
    idx = wod.profile == stations(f);


    %% interpolate (if profile doesn't overlap with glodap profile)
    % same day and within about 1km
    overlap_idx = any(abs(glodap_data.LAT-mean(wod.lat(idx))) < 0.01 & ...
        abs(glodap_data.LON-mean(wod.lon(idx))) < 0.01 & ...
        abs(glodap_data.TIME-mean(wod.time(idx))) == 0);
    
    if ~overlap_idx

        % loop through station ids, interpolate profiles and log data
        for k = 1:numel(vars) % for each variable
    
            % get temporary pressure
            temp_pres = wod.pres(idx);
            temp_var = wod.(vars{k})(idx);
            [~,unique_idx_pres] = unique(temp_pres);
    
            if length(unique_idx_pres) > 5 % if more than five depths are available
    
                % interpolate to edges (without extrapolation)
                temp_var_i = interp1(temp_pres(unique_idx_pres),...
                    temp_var(unique_idx_pres),zi,'linear');
    
                % if there is a greater than 0.1 % change per meter
                % at the bottom of the profile, remove extrapolated values
                pres_axis = temp_pres(unique_idx_pres);
                var_axis = temp_var(unique_idx_pres);
                if abs((100*((var_axis(end)-var_axis(end-1))./var_axis(end)))./...
                        (pres_axis(end)-pres_axis(end-1))) > 0.5
                    temp_var_i(zi>max(pres_axis)) = NaN;
                end
    
                % log interpolated profiles in data structure
                wod_data.(varsi{k}) = [wod_data.(varsi{k});temp_var_i];
    
                % clean up
                clear temp_var_i
    
            else
    
                wod_data.(varsi{k}) = [wod_data.(varsi{k});nan(length(zi),1)];
    
            end
    
            % clean up
            clear k temp_var temp_pres unique_idx_pres
    
        end

    %% log extra data in interpolated data structure
    wod_data.LAT = [wod_data.LAT;...
        repmat(mean(wod.lat(idx)),length(zi),1)];
    wod_data.LON = [wod_data.LON;...
        repmat(mean(wod.lon(idx)),length(zi),1)];
    wod_data.PRES = [wod_data.PRES;zi];
    wod_data.TIME = [wod_data.TIME;...
        repmat(mean(wod.time(idx)),length(zi),1)];
    wod_data.YEAR = [wod_data.YEAR;...
        repmat(mean(wod.year(idx)),length(zi),1)];
    wod_data.DAY = [wod_data.DAY;...
        repmat(mean(wod.day(idx)),length(zi),1)];
    wod_data.CRU = [wod_data.CRU;...
        repmat(mean(wod.cruise(idx)),length(zi),1)];
    wod_data.ID = [wod_data.ID;...
        repmat(mean(wod.profile(idx)),length(zi),1)];

    end

end

%% plot processed data
figure(1); hold on;
lon_temp = convert_lon(convert_lon(wod_data.LON),'0-360');
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,wod_data.LAT,'.g');
if ~exist([pwd '/' param_props.dir_name '/Figures/Data'],'dir');
    mkdir([param_props.dir_name '/Figures/Data']); end
exportgraphics(gcf,[param_props.dir_name '/Figures/Data/processed_wod_osd_filtered_' year '.png']);
close
clear lon_temp

%% Calculate absolute salinity, conservative temperature, potential density, and spice
wod_data.ABSSAL = gsw_SA_from_SP(wod_data.SAL,wod_data.PRES,wod_data.LON,wod_data.LAT);
wod_data.CNSTEMP = gsw_CT_from_t(wod_data.ABSSAL,wod_data.TEMP,wod_data.PRES);
wod_data.SIGMA = gsw_sigma0(wod_data.ABSSAL,wod_data.CNSTEMP);
wod_data.SPICE = gsw_spiciness0(wod_data.ABSSAL,wod_data.CNSTEMP);

%% display the number of casts and profiles
disp(['# of matching OSD Casts (' param_props.temp_name '): ' num2str(length(unique(wod_data.ID)))]);
disp(['# of matching Cruises (' param_props.temp_name '): ' num2str(length(unique(wod_data.CRU)))]);

%% save processed osd data
if ~exist([pwd '/' param_props.dir_name '/Data'],'dir')
    mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/processed_wod_osd_' param_props.file_name '_data_' year '.mat'],...
    'wod_data','-v7.3');

% else
% 
% % display information
% disp('WOD data already processed.')
% 
% end

end
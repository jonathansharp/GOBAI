% acquire_glodap_data
%
% DESCRIPTION:
% This function is used to import an annually updated GLODAP data file,
% then filter the data by quality flag, interpolate the data onto standard
% depth levels, and save and plot the processed glodap data data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 2/5/2025

function acquire_glodap_data(param_props,glodap_year)

%% change temporary param name for DIC
if strcmp(param_props.temp_name,'PH')
    param_props.temp_name = 'DIC';
end

%% Only do all this if downloaded glodap matlab file does not exist
if exist([param_props.dir_name '/Data/processed_glodap_' param_props.file_name '_data_' num2str(glodap_year) '.mat'],'file') ~= 2

%% load GLODAP data
year = num2str(glodap_year);
glodap_data = load(['GLODAP/GLODAPv2.' year '/GLODAPv2.' year '_Merged_Master_File.mat']);
glodap_data.time = datenum([glodap_data.G2year glodap_data.G2month ...
                       glodap_data.G2day]);
glodap_data.date = datevec(glodap_data.time);
glodap_data.date0 = glodap_data.date;
glodap_data.date0(:,2:3) = 1;
glodap_data.day = datenum(glodap_data.date) - datenum(glodap_data.date0) + 1;

%% plot unprocessed profile locations
figure(1); hold on;
set(gcf,'visible','on','position',[100 100 1600 800]);
m_proj('robinson','lon',[20 380]);
m_coast('patch',rgb('gray'));
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',-90:30:90);
lon_temp = convert_lon(convert_lon(glodap_data.G2longitude));
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,glodap_data.G2latitude,'.k');
hold off;
clear lon_temp

%% indices for glodap
idx_nans = ~isnan(glodap_data.G2temperature) & ~isnan(glodap_data.G2pressure) & ...
          ~isnan(glodap_data.G2salinity) & ~isnan(glodap_data.(param_props.glodap_name));
idx_qc = glodap_data.G2salinityqc == 1 & ...
    glodap_data.([param_props.glodap_name 'qc']) == 1;
% check for good flags
idx_flags = glodap_data.G2salinityf == 2 & glodap_data.([param_props.glodap_name 'f']) == 2;
% check depth and/or time range
idx_lims = glodap_data.G2pressure <= 2500;
%idx_lims = glodap_data.G2pressure <= 2500 & glodap_data.time > datenum(2004,1,0);
% combine indices
idx = idx_nans & idx_qc & idx_flags & idx_lims;
clear idx_nans idx_qc idx_flags idx_lims

%% remove extraneous data points
if strcmp(param_props.glodap_name,'G2tco2')
    glodap_data.G2phtsinsitutp = glodap_data.G2phtsinsitutp(idx);
    glodap_data.G2talk = glodap_data.G2talk(idx);
    glodap_data.G2oxygen = glodap_data.G2oxygen(idx);
    glodap_data.G2nitrate = glodap_data.G2nitrate(idx);
elseif strcmp(param_props.glodap_name,'G2nitrate')
    glodap_data.G2oxygen = glodap_data.G2oxygen(idx);
end
glodap_data.(param_props.glodap_name) = glodap_data.(param_props.glodap_name)(idx);
glodap_data.G2latitude = glodap_data.G2latitude(idx);
glodap_data.G2longitude = glodap_data.G2longitude(idx);
glodap_data.G2pressure = glodap_data.G2pressure(idx);
glodap_data.G2temperature = glodap_data.G2temperature(idx);
glodap_data.G2salinity = glodap_data.G2salinity(idx);
glodap_data.time = glodap_data.time(idx);
glodap_data.G2year = glodap_data.G2year(idx);
glodap_data.day = glodap_data.day(idx);
glodap_data.G2cruise = double(glodap_data.G2cruise(idx));
glodap_data.G2station = glodap_data.G2station(idx);
glodap_data.G2id = glodap_data.G2cruise.*100000+glodap_data.G2station;

%% pre-allocate glodap data structure
if strcmp(param_props.glodap_name,'G2tco2')
    glodap_data.(param_props.temp_name) = [];
    glodap_data.PH = [];
    glodap_data.TA = [];
    glodap_data.OXY = [];
    glodap_data.NIT = [];
elseif strcmp(param_props.glodap_name,'G2nitrate')
    glodap_data.(param_props.temp_name) = [];
    glodap_data.OXY = [];
else
    glodap_data.(param_props.temp_name) = [];
end
glodap_data.LAT = [];
glodap_data.LON = [];
glodap_data.PRES = [];
glodap_data.TIME = [];
glodap_data.YEAR = [];
glodap_data.DAY = [];
glodap_data.TEMP = [];
glodap_data.SAL = [];
glodap_data.CRU = [];
glodap_data.ID = [];

%% construct depth axis on which to interpolate
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
% define variables to interpolate
if strcmp(param_props.glodap_name,'G2tco2')
    vars = {'G2salinity' 'G2temperature' param_props.glodap_name 'G2phtsinsitutp' 'G2talk' 'G2oxygen' 'G2nitrate'};
    varsi = {'SAL' 'TEMP' param_props.temp_name 'PH' 'TA' 'OXY' 'NIT'};
elseif strcmp(param_props.glodap_name,'G2nitrate')
    vars = {'G2salinity' 'G2temperature' param_props.glodap_name 'G2oxygen'};
    varsi = {'SAL' 'TEMP' param_props.temp_name 'OXY'};
else
    vars = {'G2salinity' 'G2temperature' param_props.glodap_name};
    varsi = {'SAL' 'TEMP' param_props.temp_name};
end
% define station ids
stations = unique(glodap_data.G2id);

%% process glodap data
for f = 1:length(stations) % for each unique station id

    %% index according to station id
    idx = glodap_data.G2id == stations(f);

    %% interpolate
    % loop through station ids, interpolate profiles and log data
    for k = 1:numel(vars) % for each variable

        % get temporary pressure
        temp_pres = glodap_data.G2pressure(idx);
        temp_var = glodap_data.(vars{k})(idx);
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
            glodap_data.(varsi{k}) = [glodap_data.(varsi{k});temp_var_i];

            % clean up
            clear temp_var_i

        else

            glodap_data.(varsi{k}) = [glodap_data.(varsi{k});nan(length(zi),1)];

        end

        % clean up
        clear k temp_var temp_pres unique_idx_pres

    end

    % if any(glodap_data.(param_props.temp_name)(:)>100000)
    %     keyboard
    % end

    %% log extra data in interpolated data structure
    glodap_data.LAT = [glodap_data.LAT;...
        repmat(mean(glodap_data.G2latitude(idx)),length(zi),1)];
    glodap_data.LON = [glodap_data.LON;...
        repmat(mean(glodap_data.G2longitude(idx)),length(zi),1)];
    glodap_data.PRES = [glodap_data.PRES;zi];
    glodap_data.TIME = [glodap_data.TIME;...
        repmat(mean(glodap_data.time(idx)),length(zi),1)];
    glodap_data.YEAR = [glodap_data.YEAR;...
        repmat(mean(glodap_data.G2year(idx)),length(zi),1)];
    glodap_data.DAY = [glodap_data.DAY;...
        repmat(mean(glodap_data.G2day(idx)),length(zi),1)];
    glodap_data.CRU = [glodap_data.CRU;...
        repmat(mean(glodap_data.G2cruise(idx)),length(zi),1)];
    glodap_data.ID = [glodap_data.ID;...
        repmat(mean(glodap_data.G2id(idx)),length(zi),1)];

    if rem(mean(glodap_data.G2cruise(idx)),1)~=0
        keyboard
    end

end

%% clean up
clear f idx stations vars varsi zi

%% Calculate absolute salinity, conservative temperature, potential density, and spice
glodap_data.ABSSAL = gsw_SA_from_SP(glodap_data.SAL,glodap_data.PRES,glodap_data.LON,glodap_data.LAT);
glodap_data.CNSTEMP = gsw_CT_from_t(glodap_data.ABSSAL,glodap_data.TEMP,glodap_data.PRES);
glodap_data.SIGMA = gsw_sigma0(glodap_data.ABSSAL,glodap_data.CNSTEMP);
glodap_data.SPICE = gsw_spiciness0(glodap_data.ABSSAL,glodap_data.CNSTEMP);

%% display the number of matching cruises and profiles
disp(['# of matching GLODAP profiles (' param_props.temp_name '): ' num2str(length(unique(glodap_data.ID)))]);
disp(['# of matching GLODAP cruises (' param_props.temp_name '): ' num2str(length(unique(glodap_data.CRU)))]);

%% plot processed profile locations on top of unprocessed
figure(1); hold on;
lon_temp = convert_lon(convert_lon(glodap_data.LON));
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,glodap_data.LAT,'.g'); hold off;
if ~exist([param_props.dir_name '/Figures/Data'],'dir'); mkdir([param_props.dir_name '/Figures/Data']); end
export_fig(gcf,[param_props.dir_name '/Figures/Data/processed_glodap_' year '.png'],'-transparent');
clear lon_temp

%% clean up
clear path idx default_names index

%% remove unprocessed data
vars = fieldnames(glodap_data);
idx = startsWith(vars,'G2') | startsWith(vars,'expocode') | strcmp(vars,'time') | ...
    strcmp(vars,'date') | strcmp(vars,'date0') | strcmp(vars,'day');
glodap_data = rmfield(glodap_data,vars(idx));
clear idx

%% save glodap data
if ~exist([param_props.dir_name '/Data'],'dir'); mkdir([param_props.dir_name '/Data']); end
save([param_props.dir_name '/Data/processed_glodap_' param_props.file_name '_data_' year '.mat'],'glodap_data');

%% clean up
clear glodap_data
close all

% display information
disp('GLODAP data processed and saved.')

else

% display information
disp('GLODAP data already processed.')

end

end
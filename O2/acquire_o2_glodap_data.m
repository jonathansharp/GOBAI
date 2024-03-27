% acquire_o2_glodap_data
%
% DESCRIPTION:
% This function is used to import an annually updated GLODAP data file,
% then filter the data by quality flag, interpolate the data onto standard
% depth levels, and save and plot the processed glodap data data.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 09/12/2023

function acquire_o2_glodap_data(glodap_year)

%% Only do all this if downloaded glodap matlab file does not exist
% if exist(['Data/processed_glodap_o2_data_' num2str(glodap_year) '.mat'],'file') ~= 2

%% load GLODAP data
year = num2str(glodap_year);
glodap_data = load(['GLODAP/GLODAPv2.' year '/GLODAPv2.' year '_Merged_Master_File.mat']);
glodap_data.time = datenum([glodap_data.G2year glodap_data.G2month ...
                       glodap_data.G2day]);
glodap_data.date = datevec(glodap_data.time);
glodap_data.date0 = glodap_data.date;
glodap_data.date0(:,2:3) = 1;
glodap_data.day = datenum(glodap_data.date) - datenum(glodap_data.date0) + 1;

%% plot unprocessed oxygen profile locations
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

%% indices for glodap oxygen
idx_nans = ~isnan(glodap_data.G2temperature) & ~isnan(glodap_data.G2pressure) & ...
          ~isnan(glodap_data.G2salinity) & ~isnan(glodap_data.G2oxygen);
idx_qc = glodap_data.G2salinityqc == 1 & glodap_data.G2oxygenqc == 1;
idx_flags = glodap_data.G2salinityf == 2 & glodap_data.G2oxygenf == 2;
idx_lims = glodap_data.G2pressure <= 2500 & glodap_data.time > datenum(2004,1,0);
idx_oxy = idx_nans & idx_qc & idx_flags & idx_lims;
clear idx_nans idx_qc idx_flags idx_lims

%% remove extraneous data points
glodap_data.G2oxygen = glodap_data.G2oxygen(idx_oxy);
glodap_data.G2latitude = glodap_data.G2latitude(idx_oxy);
glodap_data.G2longitude = glodap_data.G2longitude(idx_oxy);
glodap_data.G2pressure = glodap_data.G2pressure(idx_oxy);
glodap_data.G2temperature = glodap_data.G2temperature(idx_oxy);
glodap_data.G2salinity = glodap_data.G2salinity(idx_oxy);
glodap_data.time = glodap_data.time(idx_oxy);
glodap_data.G2year = glodap_data.G2year(idx_oxy);
glodap_data.day = glodap_data.day(idx_oxy);
glodap_data.G2cruise = double(glodap_data.G2cruise(idx_oxy));
glodap_data.G2station = glodap_data.G2station(idx_oxy);
glodap_data.G2id = glodap_data.G2cruise.*100000+glodap_data.G2station;

%% pre-allocate glodap data structure
glodap_data.OXY = [];
glodap_data.OXY_LAT = [];
glodap_data.OXY_LON = [];
glodap_data.OXY_PRES = [];
glodap_data.OXY_TIME = [];
glodap_data.OXY_YEAR = [];
glodap_data.OXY_DAY = [];
glodap_data.OXY_TEMP = [];
glodap_data.OXY_SAL = [];
glodap_data.OXY_CRU = [];
glodap_data.OXY_ID = [];

%% construct depth axis on which to interpolate
zi = ([2.5 10:10:170 182.5 200:20:440 462.5 500:50:1350 1412.5 1500:100:1900 1975])';
% define variables to interpolate
vars = {'G2salinity' 'G2temperature' 'G2oxygen'};
varsi = {'OXY_SAL' 'OXY_TEMP' 'OXY'};
% define station ids
stations = unique(glodap_data.G2id);

%% process glodap oxygen data
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

    if any(glodap_data.OXY(:)>1000)
        keyboard
    end

    %% log extra data in interpolated data structure
    glodap_data.OXY_LAT = [glodap_data.OXY_LAT;...
        repmat(mean(glodap_data.G2latitude(idx)),length(zi),1)];
    glodap_data.OXY_LON = [glodap_data.OXY_LON;...
        repmat(mean(glodap_data.G2longitude(idx)),length(zi),1)];
    glodap_data.OXY_PRES = [glodap_data.OXY_PRES;zi];
    glodap_data.OXY_TIME = [glodap_data.OXY_TIME;...
        repmat(mean(glodap_data.time(idx)),length(zi),1)];
    glodap_data.OXY_YEAR = [glodap_data.OXY_YEAR;...
        repmat(mean(glodap_data.G2year(idx)),length(zi),1)];
    glodap_data.OXY_DAY = [glodap_data.OXY_DAY;...
        repmat(mean(glodap_data.G2day(idx)),length(zi),1)];
    glodap_data.OXY_CRU = [glodap_data.OXY_CRU;...
        repmat(mean(glodap_data.G2cruise(idx)),length(zi),1)];
    glodap_data.OXY_ID = [glodap_data.OXY_ID;...
        repmat(mean(glodap_data.G2id(idx)),length(zi),1)];

    if rem(mean(glodap_data.G2cruise(idx)),1)~=0
        keyboard
    end

end

%% clean up
clear f idx stations vars varsi zi

%% Calculate absolute salinity, conservative temperature, potential density, and spice
glodap_data.OXY_ABSSAL = gsw_SA_from_SP(glodap_data.OXY_SAL,glodap_data.OXY_PRES,glodap_data.OXY_LON,glodap_data.OXY_LAT);
glodap_data.OXY_CNSTEMP = gsw_CT_from_t(glodap_data.OXY_ABSSAL,glodap_data.OXY_TEMP,glodap_data.OXY_PRES);
glodap_data.OXY_SIGMA = gsw_sigma0(glodap_data.OXY_ABSSAL,glodap_data.OXY_CNSTEMP);
glodap_data.OXY_SPICE = gsw_spiciness0(glodap_data.OXY_ABSSAL,glodap_data.OXY_CNSTEMP);

%% display the number of matching cruises and profiles
disp(['# of matching GLODAP profiles (OXY): ' num2str(length(unique(glodap_data.OXY_ID)))]);
disp(['# of matching GLODAP cruises (OXY): ' num2str(length(unique(glodap_data.OXY_CRU)))]);

%% plot processed oxygen profile locations on top of unprocessed
figure(1); hold on;
lon_temp = convert_lon(convert_lon(glodap_data.OXY_LON));
lon_temp(lon_temp < 20) = lon_temp(lon_temp < 20) + 360;
m_scatter(lon_temp,glodap_data.OXY_LAT,'.g'); hold off;
if ~exist('Figures','dir'); mkdir('Figures'); end
exportgraphics(gcf,['O2/Figures/Data/processed_glodap_' year '.png']);
clear lon_temp

%% clean up
clear path idx_oxy default_names index

%% remove unprocessed data
vars = fieldnames(glodap_data);
idx = startsWith(vars,'OXY');
glodap_data = rmfield(glodap_data,vars(~idx));
clear idx

%% save glodap oxygen data
if ~exist([pwd '/O2/Data'],'dir'); mkdir('/O2/Data'); end
save(['O2/Data/processed_glodap_o2_data_' year '.mat'],'glodap_data');

%% clean up
clear glodap_data
close all

% end

end
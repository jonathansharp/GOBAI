% define mask for gobai files
% ver: 'v1.0-HR', 'v1.1-HR'
% var: 'O2', 'NO3', 'DIC'

function mask = gobai_mask(ver,ext,var)

if strcmp(var,'O2')
    varname = 'o2';
elseif strcmp(var,'NO3')
    varname = 'no3';
elseif strcmp(var,'DIC')
    varname = 'dic';
end

 % file path
path = ['/raid/Data/GOBAI-' var '/' ver '/'];

% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'longitude');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'latitude');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'mean_pressure');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver ext '.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);

% process time axis
if strcmp(ver,'v1.0') || strcmp(ver,'v2.0')
    time_axis = double(GOBAI.time);
elseif strcmp(ver,'v2.1')
    time_axis = datenum(2004,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v2.2') || strcmp(ver,'v2.3')
    time_axis = datenum(1950,0,0) + double(GOBAI.time);
elseif strcmp(ver,'v1.0-HR') || strcmp(ver,'v1.1-HR')
    time_axis = datenum(1950,1,1) + double(GOBAI.time);
end

if strcmp(ver,'v1.0-HR') || strcmp(ver,'v1.1-HR')
    GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
    if ~isfile([path 'mask-' ver ext '.mat'])
        mask = true(size(GOBAI.vol));
        for t = 1:length(GOBAI.time)
            gobai_tmp = ncread([path 'GOBAI-' var '-' ver ext '.nc'],...
                varname,[1 1 1 t],[Inf Inf Inf 1]);
            try
            mask(isnan(gobai_tmp)) = false;
            catch
            keyboard
            end
        end
        save([path 'mask-' ver ext '.mat'],'mask');
    else
        load([path 'mask-' ver ext '.mat'],'mask');
    end
end

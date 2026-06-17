% calculate anomalies for 
% ver: 'v1.0-HR', 'v1.1-HR'
% var: 'O2', 'NO3', 'DIC'

function mask = gobai_fit_models(ver,var)

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
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');

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

% fit models
coeffs.x = cell(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres));
for x = 1%:length(GOBAI.lon)
    for y = 1:length(GOBAI.lat)
        for z = 1%:length(GOBAI.pres)
            gobai_tmp = squeeze(ncread([path 'GOBAI-' var '-' ver '.nc'],...
                 varname,[x y z 1],[1 1 1 Inf]));
            if sum(~isnan(gobai_tmp)) >= 500
                [~,~,coeffs.x{x,y,z}] = leastsq2(time_axis,gobai_tmp,time_axis(1),2,[365.2425 365.2425/2])
            else
            
            end
        end
    end
end

mask(isnan(gobai_tmp)) = false;

save([path 'mask-' ver '.mat'],'mask');

load([path 'mask-' ver '.mat'],'mask');


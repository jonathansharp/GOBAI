%% v 1.0-HR-FFNN
% file information
ver = 'v1.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '-HR/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],'time');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
% define mask
if ~isfile([path 'mask-' ver '.mat'])
    mask = true(size(GOBAI.vol));
    for t = 1:length(GOBAI.time)
        gobai_tmp = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],...
            'o2',[1 1 1 t],[Inf Inf Inf 1]);
        mask(isnan(gobai_tmp)) = false;
    end
    save([path 'mask-' ver '.mat'],'mask');
else
    load([path 'mask-' ver '.mat']);
end
GOBAI.vol(~mask) = NaN;
% download oxygen and calculate global mean and inventory
global_mean = nan(size(GOBAI.time));
global_inv = nan(size(GOBAI.time));
for t = 1:length(GOBAI.time)
    gobai_tmp = ncread([path 'GOBAI-' var '-' ver '-HR.nc'],...
        'o2',[1 1 1 t],[Inf Inf Inf 1]);
    gobai_tmp(~mask) = NaN;
    % umol/kg * kg/m3 * m3
    global_inv(t) = sum(gobai_tmp(:).*GOBAI.vol(:)
    global_mean(t) = sum(gobai_tmp(:).*...
        GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
end
% calculate global mean
plot(datenum(1950,0,0)+double(GOBAI.time),global_mean,'LineWidth',2);
save([path 'mask'],'mask');
save([path 'global_mean'],'global_mean');
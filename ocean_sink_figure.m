% GCB Sink Figure

% gcb data
gcb_f = 'Data/Global_Carbon_Budget_2025_v1.0.xlsx'; % filename from 2025 budget
os = readtable(gcb_f,'Sheet','Ocean Sink','NumHeaderLines',31); % fluxes in PgC

% gobai data
ver = 'HR-v1.1'; % product version
path = ['/raid/Data/GOBAI-DIC/' ver '/']; % file path
gobai_f = [path 'GOBAI-Monthly-Mean-1x1-DIC-HR-v202602.nc']; % file name
rfrom_t_f = '/raid/Data/RFROM/RFROM_TEMP_v2.2-2025/RFROMV22_TEMP_MONTHLY_MEAN_1x1.nc';
rfrom_s_f = '/raid/Data/RFROM/RFROM_SAL_v2.2-2025/RFROMV22_SAL_MONTHLY_MEAN_1x1.nc';
% download dimensions
gobai.lon = ncread(gobai_f,'longitude');
gobai.lat = ncread(gobai_f,'latitude');
gobai.pres = ncread(gobai_f,'mean_pressure');
gobai.time = datenum(1950,1,1) + ncread(gobai_f,'time');
% calculate mass per grid cell
gobai.vol = weights3d(gobai.lon,gobai.lat,gobai.pres);
% determine mask (time varying coverage)
mask = true(size(gobai.vol));
for t = 1:length(gobai.time)
    dic = ncread(gobai_f,'dic',[1 1 1 t],[Inf Inf Inf 1]);
    mask(isnan(dic)) = false;
end

% determine inventory (PgC)
[lon3d,lat3d,pres3d] = ndgrid(gobai.lon,gobai.lat,gobai.pres);
year = datevec(double(gobai.time)); year = year(:,1);
dic_inv = nan(length(gobai.time),1);
canth_inv = nan(length(gobai.time),1);
canth_3d = nan(size(lon3d));
for t = 1:length(gobai.time)
    dic = ncread(gobai_f,'dic',[1 1 1 t],[Inf Inf Inf 1]);
    tmp = ncread(rfrom_t_f,'ocean_temperature',[1 1 1 t],[Inf Inf Inf 1]);
    sal = ncread(rfrom_s_f,'ocean_salinity',[1 1 1 t],[Inf Inf Inf 1]);
    dens = gsw_rho(sal,tmp,pres3d); % kg/m3
    % calculate global dic inventory 
    dic_inv(t) = ... % Pg = umol/kg * kg/m3 * m3 * g/mol * mol/umol
        sum(dic(:).*dens(:).*gobai.vol(:).*12.01.*10^-21,'omitnan');
    % apply TRACE to calculate canth
    lon2d = lon3d(:,:,1); lat2d = lat3d(:,:,1);
    for p = 1:length(gobai.pres)
        pres2d = pres3d(:,:,p); tmp2d = tmp(:,:,p); sal2d = sal(:,:,p);
        canth_2d = nan(size(sal2d));
        canth_2d(mask(:,:,p)) = ...
            TRACEv1([lon2d(mask(:,:,p)) lat2d(mask(:,:,p)) ...
            pres2d(mask(:,:,p))],year(t),[tmp2d(mask(:,:,p)) ...
            sal2d(mask(:,:,p))],[2 1],1);
        canth_3d(:,:,p) = canth_2d;
    end
    % calculate global canth inventory 
    canth_inv(t) = ... % Pg = umol/kg * kg/m3 * m3 * g/mol * mol/umol
        sum(canth_3d(:).*dens(:).*gobai.vol(:).*12.01.*10^-21,'omitnan');
    % print status
    % gobai.time
end

% annual mean of dic inventory
dic_inv_ann = mean(reshape(gobai.dic_inv,12,[]),1);
time_ann = mean(reshape(gobai.time,12,[]),1);

% scaling factor
scal = 38000./mean(gobai.dic_inv);

% plot ocean carbon sink figure
figure; plot(os.year,os.GCB);
figure; plot(gobai.time,gobai.dic_inv.*scal);
figure; plot(double(gobai.time),movmean(dic_inv,12)); datetick('x','yyyy')
figure; plot(double(time_ann(2:end)),diff(dic_inv_ann.*scal)); datetick('x','yyyy')
figure; plot(double(time_ann),dic_inv_ann.*scal); datetick('x','yyyy')

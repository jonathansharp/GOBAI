% file name
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc';
fname = '/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc';

% download dimensions
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
pres = ncread(fname,'pres');
time = ncread(fname,'time');

% [vol,area,h] = weights3d(lon,lat,pres);

% get long-term mean (1993-2023)
gobai_mean_lt = nan(length(lon),length(lat));
for t = 1:length(time)
    date = datevec(double(time(t)+datenum(1950,0,1))); year = date(1);
    if year < 2024
        gobai_tmp = ncread(fname,'oxy',[1 1 1 t],[Inf Inf Inf 1]);
        gobai_mean_lt = moving_mean(gobai_mean_lt,gobai_tmp,t);
        if date(2) == 1 && date(3) <= 7
            disp(['Starting ' num2str(year)]);
        end
    end
end

% get 2024 mean
years = datevec(double(time+datenum(1950,0,1))); years = years(:,1);
year_idx = find(years == 2024);
gobai_tmp = ncread(fname,'oxy',[1 1 1 year_idx(1)],[Inf Inf Inf length(year_idx)]);
gobai_mean_2024 = mean(gobai_tmp,4,'omitnan');

% calculate inventories
[~,~,height] = weights3d(lon,lat,pres);
idx_both = ~isnan(gobai_mean_lt) & ~isnan(gobai_mean_2024);
height(~idx_both) = NaN; gobai_mean_lt(~idx_both) = NaN; gobai_mean_2024(~idx_both) = NaN;
gobai_inv_lt = sum(gobai_mean_lt.*height,3,'omitnan')./(10.^3); % mmol/m2
gobai_inv_2024 = sum(gobai_mean_2024.*height,3,'omitnan')./(10.^3); % mmol/m2

% plot anomaly
figure('position',[100 100 800 500]);
worldmap([-70 70],[20 380]);
pcolorm(double(lat),double(lon),(gobai_inv_2024 - gobai_inv_lt)');
colorbar(gca,'southoutside');
clim([-15 15]);
colormap(cmocean('balance','pivot',0));
plot_land('map');
title('Oxygen content anomaly from 0 to 2000 dbar in 2024 (mmol m^{-2})')
plabel off; mlabel off;
export_fig(gcf,'Figures/o2_anom_v1.0-HR.png','-transparent');
close


clear
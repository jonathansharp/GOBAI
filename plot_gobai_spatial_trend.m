%fname = '/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc';
%fname = '/raid/Data/GOBAI-O2/v2.2/GOBAI-O2-v2.2.nc';
%fname = '/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc';
%fname = '/raid/Data/GOBAI-O2/v2.4/GOBAI-O2-v2.4-two-depth-adjustment.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc';
fname = '/raid/Data/GOBAI-O2/HR-v1.3/GOBAI-Monthly-Mean-1x1-O2-HR-v202606.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-FFNN-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-GBM-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-RFR-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-NO3/v1.1-HR/GOBAI-NO3-v1.1-HR.nc';
%fname = '/raid/Data/GOBAI-DIC/v1.0-HR/GOBAI-DIC-FFNN-v1.0-HR.nc';

lon = ncread(fname,'longitude');
lat = ncread(fname,'latitude');
pres = ncread(fname,'mean_pressure');
time = ncread(fname,'time');
%gobai = squeeze(ncread(fname,'oxy',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'oxygen',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = ncread(fname,'o2');
%gobai = squeeze(ncread(fname,'no3',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'dic',[1 1 11 1],[Inf Inf 1 Inf]));

[vol,area,h] = weights3d(lon,lat,pres);

tr = nan(length(lon),length(lat));
for x = 1:length(lon)
    for y = 1:length(lat)
        gobai_temp = squeeze(ncread(fname,'o2',[x,y,11,1],[1 1 26 Inf]));
        wts = squeeze(vol(x,y,11:36));
        col_mean = sum(gobai_temp(~isnan(gobai_temp)).*wts(~isnan(gobai_temp)))./...
            sum(wts(~isnan(gobai_temp)));
        mdl = polyfit(time,squeeze(col_mean(x,y,:)),1);
        tr(x,y) = mdl(1)*365;
    end
end

figure('position',[100 100 800 500]);
worldmap([-90 90],[20 380]);
pcolorm(double(lat),double(lon),tr');
colorbar(gca,'southoutside');
clim([-1 1]);
colormap(cmocean('balance','pivot',0));
plot_land('map');
plabel off; mlabel off;
export_fig(gcf,'Figures/o2_trend_HR_v1.3.png');
close


clear
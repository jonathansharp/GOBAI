%fname = '/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc';
%fname = '/raid/Data/GOBAI-O2/v2.2/GOBAI-O2-v2.2.nc';
fname = '/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-FFNN-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-GBM-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-O2/v1.0-HR/GOBAI-O2-RFR-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-NO3/v1.0-HR/GOBAI-NO3-FFNN-v1.0-HR.nc';
%fname = '/raid/Data/GOBAI-DIC/v1.0-HR/GOBAI-DIC-FFNN-v1.0-HR.nc';

lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
pres = ncread(fname,'pres');
time = ncread(fname,'time');
gobai = squeeze(ncread(fname,'oxy',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'oxygen',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'o2',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'no3',[1 1 11 1],[Inf Inf 1 Inf]));
%gobai = squeeze(ncread(fname,'dic',[1 1 11 1],[Inf Inf 1 Inf]));

% [vol,area,h] = weights3d(lon,lat,pres);

tr = nan(length(lon),length(lat));
for x = 1:length(lon)
    for y = 1:length(lat)
        % gobai_temp = squeeze(ncread(fname,'oxygen',[x,y,11,1],[1 1 1 Inf]));
        % wts = squeeze(vol(x,y,:));
        % col_mean = sum(gobai_temp(~isnan(gobai_temp)).*wts(~isnan(gobai_temp)))./...
        %     sum(wts(~isnan(gobai_temp)));
        mdl = polyfit(time,squeeze(gobai(x,y,:)),1);
        tr(x,y) = mdl(1)*365;
    end
end

figure('position',[100 100 800 500]);
worldmap([-90 90],[20 380]);
pcolorm(double(lat),double(lon),tr');
colorbar(gca,'southoutside');
clim([-2 2]);
colormap(cmocean('balance','pivot',0));
plot_land('map');
plabel off; mlabel off;
export_fig(gcf,'Figures/o2_trend_v2.2.png','-transparent');
close


clear
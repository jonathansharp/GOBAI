%fname1 = '/raid/Data/GOBAI-O2/v2.1/GOBAI-O2-v2.1.nc';
fname1 = '/raid/Data/GOBAI-O2/v2.2/GOBAI-O2-v2.2.nc';
%fname1 = '/raid/Data/GOBAI-O2/v3.0/GOBAI-O2-v3.0.nc';
fname2 = '/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc';


lon = ncread(fname1,'lon');
lat = ncread(fname1,'lat');
pres = ncread(fname1,'pres');
time1 = ncread(fname1,'time');
time2 = ncread(fname2,'time');
gobai1 = squeeze(ncread(fname1,'uncer',[1 1 11 1],[Inf Inf 1 Inf]));
gobai2= squeeze(ncread(fname2,'uncer',[1 1 11 1],[Inf Inf 1 Inf]));

gobai1 = squeeze(mean(mean(gobai1,1,'omitnan'),2,'omitnan'));
gobai2 = squeeze(mean(mean(gobai2,1,'omitnan'),2,'omitnan'));

figure('position',[100 100 800 500]);
worldmap([-90 90],[20 380]);
pcolorm(double(lat),double(lon),mean(gobai1,3,'omitnan')');
pcolorm(double(lat),double(lon),mean(gobai1,3,'omitnan')'-mean(gobai2,3,'omitnan')');
colorbar(gca,'southoutside');
clim([-20 20]);
colormap(cmocean('balance','pivot',0));
plot_land('map');
plabel off; mlabel off;
export_fig(gcf,'Figures/no3_trend.png','-transparent');
%close


clear
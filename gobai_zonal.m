%% Plot global zonal mean [O2]

% initialize figure
figure;
set(gcf,'units','inches','position',[0 14 20 10]); hold on;
set(gca,'fontsize',26,'YDir','reverse','linewidth',3,'layer','top','box','on');
title('GOBAI-O_{2} Zonal Mean [O_{2}] (\mumol kg^{-1})','fontsize',32);

% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate 2d mean
gobai_oxy_3d = mean(GOBAI.oxy,4,'omitnan');
global_mean_2d = sum(gobai_oxy_3d.*GOBAI.vol,1,'omitnan')./sum(GOBAI.vol,1,'omitnan');
global_mean_2d = squeeze(global_mean_2d);

% plot data
contourf(GOBAI.lat,GOBAI.pres,global_mean_2d',38.5:25:362.5);

% modify axes
xlim([-64.5 64.5]);
% ylim([0 1500]);
xlabel('Latitude');
ylabel('Pressure (dbars)');
caxis([38.5 362.5]);

% colorbar
c=colorbar;
c.TickLength = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
colormap(cmocean('ice',13));

% save figure
exportgraphics(gcf,['Figures/zonal_mean_O2_GOBAI_' ver '.png'])
close
clear

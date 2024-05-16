%% Plot global column inventory O2

% initialize figure
figure; worldmap([-64.5 79.5],[20 380]); hold on;
set(gcf,'units','inches','position',[0 14 20 10]);
framem('FlineWidth',3,'FEdgeColor','k');
setm(gca,'fontsize',26,'ffacecolor','w');
set(gca,'layer','top');
title('GOBAI-O_{2} Column Mean [O_{2}] (\mumol kg^{-1})','fontsize',32);

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
global_mean_2d = sum(gobai_oxy_3d.*GOBAI.vol,3,'omitnan')./sum(GOBAI.vol,3,'omitnan');

% plot data
contourfm(double(GOBAI.lat),double(GOBAI.lon),double(global_mean_2d)',-12.5:25:362.5);
land = shaperead('landareas', 'UseGeoCoords',true);
geoshow(land,'FaceColor',rgb('grey'));
clim([-12.5 362.5]);
mlabel off; plabel off;

% colorbar
c=colorbar;
c.TickLength  = 0;
c.LineWidth = 3;
c.Position(3) = 2*c.Position(3);
c.Position(1) = c.Position(1)+0.04;
c.FontSize = 26;
colormap(cmocean('ice',17));

% save figure
exportgraphics(gcf,['Figures/col_mean_O2_GOBAI_' ver '.png']);
close
clear

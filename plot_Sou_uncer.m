%% load GOBAI v2.2
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
GOBAI.oxy_uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');

%% depth level
d = 120; index = find(GOBAI.pres == d);

%% plot
figure; worldmap([-90 -30],[20 380]); hold on;
set(gcf,'units','normalized','position',[0 1 0.5 0.5]);
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
pcolorm(double(GOBAI.lat),double([GOBAI.lon;GOBAI.lon(end)+1]),...
    double([mean(GOBAI.oxy_uncer(:,:,index,:),4);...
    mean(GOBAI.oxy_uncer(end,:,index,:),4)])');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
c=colorbar; caxis([8 20]);
colormap(cmocean('deep',12));
c.FontSize = 14; c.TickLength = 0;
c.Label.String = '[O_{2}] Uncer. [\mumol kg^{-1}]';
mlabel off; plabel off;
exportgraphics(gcf,['Figures/Sou_ESM4_oxy_uncer_' num2str(d) '.png']);
close

%% clean up
clear
